# gis_visualizer_41.py
import math
import os
import tempfile
import time
from functools import lru_cache

import bpy
import geopandas as gpd
import numpy as np
import rasterio
from bpy.props import (
    BoolProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    PointerProperty,
    StringProperty,
)
from bpy.types import EnumProperty, Operator, Panel, PropertyGroup
from mathutils import Matrix, Vector
from rasterio import features
from shapely.geometry import MultiPolygon, Polygon

bl_info = {
    "name": "GIS Visualizer",
    "author": "GIS Team",
    "version": (1, 0),
    "blender": (4, 1, 0),
    "location": "View3D > Sidebar > 工具",
    "description": "导入并可视化GIS数据",
    "category": "3D View",
}


class GisProps(PropertyGroup):
    # 矢量数据属性
    shp_path: StringProperty(
        name="矢量文件", subtype="FILE_PATH", description="选择SHP矢量文件"
    )
    height_field: StringProperty(
        name="高度字段", default="height", description="用于确定挤出高度的属性字段"
    )
    name_field: StringProperty(
        name="名称字段", default="name", description="用于命名对象的属性字段"
    )
    height_scale: FloatProperty(
        name="高度缩放",
        default=1.0,
        min=0.01,
        max=100.0,
        description="控制挤出高度的缩放因子",
    )
    height_offset: FloatProperty(
        name="高度基准", default=0.0, description="Z轴基准面偏移值"
    )
    auto_center: BoolProperty(
        name="自动居中", default=True, description="自动将数据中心点置于场景原点"
    )

    # 栅格数据属性
    tif_path: StringProperty(
        name="栅格文件", subtype="FILE_PATH", description="选择TIF栅格文件"
    )
    raster_scale: FloatProperty(
        name="栅格高度缩放",
        default=1.0,
        min=0.01,
        max=100.0,
        description="控制栅格高度的缩放因子",
    )
    # cube_size: FloatProperty(
    #     name="栅格尺寸",
    #     default=0.9,
    #     min=0.1,
    #     max=1.0,
    #     description="控制每个栅格单元的尺寸",
    # )
    gap_size: FloatProperty(
        name="间隙比例",
        default=0.05,
        min=0.0,
        max=0.3,
        description="控制栅格单元之间的间隙大小",
    )

    # # 内存优化选项
    # chunk_size: IntProperty(
    #     name="分块大小",
    #     default=256,
    #     min=64,
    #     max=1024,
    #     description="栅格数据处理的分块大小",
    # )
    use_lod: BoolProperty(
        name="使用LOD", default=True, description="启用详细程度（LOD）系统来优化性能"
    )
    lod_distance: FloatProperty(
        name="LOD距离", default=100.0, min=10.0, max=1000.0, description="LOD切换距离"
    )


class GIS_OT_ImportVector(Operator):
    bl_idname = "gis.import_vector"
    bl_label = "导入矢量数据"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        props = context.scene.gis_props
        if not props.shp_path or not os.path.exists(props.shp_path):
            self.report({"ERROR"}, "请选择有效的SHP文件")
            return {"CANCELLED"}

        try:
            start_time = time.time()
            self.process_vector(context, props)
            self.report(
                {"INFO"}, f"矢量数据导入完成，耗时: {time.time() - start_time:.2f}秒"
            )
            return {"FINISHED"}
        except Exception as e:
            self.report({"ERROR"}, f"矢量处理失败: {str(e)}")
            return {"CANCELLED"}

    def get_field_names(self, gdf):
        """获取并返回GeoDataFrame的所有字段名称"""
        return list(gdf.columns)

    def process_vector(self, context, props):
        # 创建一个新的集合来存放矢量数据
        collection_name = os.path.basename(props.shp_path).split(".")[0]
        if collection_name in bpy.data.collections:
            vector_collection = bpy.data.collections[collection_name]
        else:
            vector_collection = bpy.data.collections.new(collection_name)
            context.scene.collection.children.link(vector_collection)

        # 读取SHP文件
        gdf = gpd.read_file(props.shp_path)

        # 获取坐标边界用于居中
        if props.auto_center:
            bounds = gdf.total_bounds
            center_x = (bounds[0] + bounds[2]) / 2
            center_y = (bounds[1] + bounds[3]) / 2
            offset = Vector((-center_x, -center_y, 0))
        else:
            offset = Vector((0, 0, 0))

        # 逐个处理几何要素
        for idx, row in gdf.iterrows():
            geom = row.geometry
            if geom is None:
                continue

            # 处理不同类型的几何体
            if geom.geom_type == "Polygon":
                self.create_polygon(vector_collection, geom, row, idx, props, offset)
            elif geom.geom_type == "MultiPolygon":
                for i, poly in enumerate(geom.geoms):
                    self.create_polygon(
                        vector_collection, poly, row, f"{idx}_{i}", props, offset
                    )
            # 其他类型后续可添加支持

    def create_polygon(self, collection, geom, row, idx, props, offset):
        # 获取名称
        if props.name_field and props.name_field in row:
            object_name = str(row[props.name_field])
        else:
            object_name = f"poly_{idx}"

        # 处理重名
        if object_name in bpy.data.objects:
            object_name = f"{object_name}_{idx}"

        # 获取高度值
        height = 0
        if props.height_field and props.height_field in row:
            try:
                height = float(row[props.height_field]) * props.height_scale
            except (ValueError, TypeError):
                height = 0

        # 创建顶点
        verts = []
        if isinstance(geom, Polygon):
            # 添加外环
            exterior_coords = list(geom.exterior.coords)
            verts.extend(
                [
                    Vector(
                        (coord[0] + offset.x, coord[1] + offset.y, props.height_offset)
                    )
                    for coord in exterior_coords
                ]
            )

            # 创建面
            faces = [
                list(range(len(exterior_coords) - 1))
            ]  # 最后一个点与第一个点重合，排除

            # 处理内环（孔洞）
            holes = []
            for interior in geom.interiors:
                start_idx = len(verts)
                interior_coords = list(interior.coords)
                verts.extend(
                    [
                        Vector(
                            (
                                coord[0] + offset.x,
                                coord[1] + offset.y,
                                props.height_offset,
                            )
                        )
                        for coord in interior_coords
                    ]
                )
                holes.append(
                    list(range(start_idx, start_idx + len(interior_coords) - 1))
                )

        # 创建网格
        mesh = bpy.data.meshes.new(name=object_name)
        mesh.from_pydata(verts, [], faces)
        mesh.validate()

        # 创建对象
        obj = bpy.data.objects.new(object_name, mesh)
        collection.objects.link(obj)

        # 挤出处理
        if height != 0:
            solidify = obj.modifiers.new("Solidify", "SOLIDIFY")
            solidify.thickness = abs(height)
            # solidfy offset is -1 up and 1 down
            solidify.offset = -1.0 if height > 0 else 1.0

        # 添加自定义属性存储原始数据
        for key, value in row.items():
            if key != "geometry" and not key.startswith("_"):
                obj["gis_" + key] = value


class GIS_OT_ImportRaster(Operator):
    bl_idname = "gis.import_raster"
    bl_label = "导入栅格数据"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        props = context.scene.gis_props
        if not props.tif_path or not os.path.exists(props.tif_path):
            self.report({"ERROR"}, "请选择有效的TIF文件")
            return {"CANCELLED"}

        try:
            start_time = time.time()
            self.process_raster(context, props)
            self.report(
                {"INFO"}, f"栅格数据导入完成，耗时: {time.time() - start_time:.2f}秒"
            )
            return {"FINISHED"}
        except Exception as e:
            self.report({"ERROR"}, f"栅格处理失败: {str(e)}")
            import traceback

            traceback.print_exc()
            return {"CANCELLED"}

    def process_raster(self, context, props):
        # 创建一个新的集合来存放栅格数据
        collection_name = os.path.basename(props.tif_path).split(".")[0]
        if collection_name in bpy.data.collections:
            raster_collection = bpy.data.collections[collection_name]
        else:
            raster_collection = bpy.data.collections.new(collection_name)
            context.scene.collection.children.link(raster_collection)

        # 读取TIF文件
        with rasterio.open(props.tif_path) as src:
            # 获取栅格信息
            data = src.read(1)
            transform = src.transform
            nodata = src.nodata if src.nodata is not None else -9999
            width = src.width
            height = src.height

            # 计算世界坐标系到本地坐标系的变换
            # 获取栅格的地理范围
            bounds = src.bounds
            center_x = (bounds.left + bounds.right) / 2
            center_y = (bounds.bottom + bounds.top) / 2

            # 如果启用了自动居中，则计算偏移量
            if props.auto_center:
                offset_x = -center_x
                offset_y = -center_y
            else:
                offset_x = 0
                offset_y = 0

            # 计算像素分辨率

            res_x, res_y = transform.a, abs(transform.e)
            cell_size = res_x * (1 - props.gap_size)
            # half_cell = cell_size * 0.5
            # 生成有效点云数据
            rows, cols = data.shape

            # 生成行列索引网格
            row_indices, col_indices = np.meshgrid(
                np.arange(rows), np.arange(cols), indexing="ij"
            )

            # 过滤掉无效值（nodata）
            valid_mask = data != nodata
            valid_rows = row_indices[valid_mask]
            valid_cols = col_indices[valid_mask]
            valid_values = data[valid_mask]

            # 计算有效点的地理坐标
            x_geo, y_geo = transform * (valid_cols + 0.5, valid_rows + 0.5)

            # 转换为Blender坐标（Y轴翻转）
            x = x_geo
            y = y_geo
            z = valid_values * props.raster_scale

            # 将点云数据存储为Numpy数组
            points = np.column_stack((x, y, np.zeros_like(z)))
            heights = z

            # 创建点云对象
            mesh = bpy.data.meshes.new("RasterPoints")
            mesh.vertices.add(len(points))
            mesh.vertices.foreach_set("co", [c for p in points for c in p])

            # 添加高度属性
            height_attr = mesh.attributes.new(
                name="height", type="FLOAT", domain="POINT"
            )
            height_attr.data.foreach_set("value", heights)

            # 创建点云对象
            point_cloud = bpy.data.objects.new("RasterPointCloud", mesh)
            raster_collection.objects.link(point_cloud)
            point_cloud.location = (offset_x, offset_y, 0)
            self.setup_geometry_nodes(
                point_cloud,
                cell_size,
                props,
            )
            # if props.chunk_size < width or props.chunk_size < height:
            #     # 使用几何节点进行分块处理
            #     self.create_chunked_geometry_nodes(
            #         context,
            #         raster_collection,
            #         data,
            #         nodata,
            #         transform,
            #         props,
            #         offset_x,
            #         offset_y,
            #         res_x,
            #         res_y,
            #         width,
            #         height,
            #         bounds,
            #     )
            # else:
            #     # 直接使用几何节点处理整个栅格
            #     self.create_geometry_nodes(
            #         context,
            #         raster_collection,
            #         data,
            #         nodata,
            #         transform,
            #         props,
            #         offset_x,
            #         offset_y,
            #         res_x,
            #         res_y,
            #         width,
            #         height,
            #         bounds,
            #     )

    # def create_chunked_geometry_nodes(
    #     self,
    #     context,
    #     collection,
    #     data,
    #     nodata,
    #     transform,
    #     props,
    #     offset_x,
    #     offset_y,
    #     res_x,
    #     res_y,
    #     width,
    #     height,
    #     bounds,
    # ):
    #     """使用分块方式创建几何节点栅格"""
    #     chunk_size = props.chunk_size

    #     # 计算需要的分块数量
    #     cols = math.ceil(width / chunk_size)
    #     rows = math.ceil(height / chunk_size)

    #     # 创建分块控制器对象
    #     controller = bpy.data.objects.new("RasterController", None)
    #     collection.objects.link(controller)

    #     # 为每个分块创建几何节点
    #     for row in range(rows):
    #         for col in range(cols):
    #             # 计算当前分块的范围
    #             x_start = col * chunk_size
    #             y_start = row * chunk_size
    #             x_end = min((col + 1) * chunk_size, width)
    #             y_end = min((row + 1) * chunk_size, height)

    #             # 提取当前分块的数据
    #             chunk_data = data[y_start:y_end, x_start:x_end]

    #             # 跳过全部为nodata的块
    #             if np.all(chunk_data == nodata):
    #                 continue

    #             # 创建分块对象
    #             chunk_name = f"Chunk_{row}_{col}"
    #             point_cloud_obj = self.create_point_cloud(
    #                 chunk_name,
    #                 chunk_data,
    #                 nodata,
    #                 res_x,
    #                 res_y,
    #                 bounds,
    #                 offset_x,
    #                 offset_y,
    #                 x_start,
    #                 y_start,
    #             )
    #             collection.objects.link(point_cloud_obj)

    #             # 设置分块对象的位置
    #             chunk_x = x_start * res_x + offset_x
    #             # chunk_y = y_start * res_y + offset_y
    #             chunk_y = bounds.top - (y_start * res_y) + offset_y
    #             point_cloud_obj.location = (chunk_x, chunk_y, 0)

    #             # 为分块创建几何节点
    #             self.setup_geometry_nodes(
    #                 point_cloud_obj,
    #                 point_cloud_obj,
    #                 props,
    #                 res_x,
    #                 res_y,
    #             )

    #             # 设置分块对象为控制器的子对象
    #             point_cloud_obj.parent = controller

    # def create_geometry_nodes(
    #     self,
    #     context,
    #     collection,
    #     data,
    #     nodata,
    #     transform,
    #     props,
    #     offset_x,
    #     offset_y,
    #     res_x,
    #     res_y,
    #     width,
    #     height,
    #     bounds,
    # ):
    #     """使用几何节点创建整个栅格"""
    #     # 创建控制器对象
    #     controller_name = "RasterController"
    #     point_cloud_obj = self.create_point_cloud(
    #         controller_name,
    #         data,
    #         nodata,
    #         res_x,
    #         res_y,
    #         bounds,
    #         offset_x,
    #         offset_y,
    #         0,
    #         0,
    #     )
    #     point_cloud_obj.location = (offset_x, offset_y, 0)
    #     collection.objects.link(point_cloud_obj)

    #     # 设置几何节点
    #     self.setup_geometry_nodes(point_cloud_obj, point_cloud_obj, props, res_x, res_y)

    # def create_point_cloud(
    #     self,
    #     obj_name,
    #     data_chunk,
    #     nodata,
    #     res_x,
    #     res_y,
    #     bounds,
    #     offset_x,
    #     offset_y,
    #     x_start,
    #     y_start,
    # ):
    #     valid_mask = data_chunk != nodata
    #     rows, cols = np.nonzero(valid_mask)

    #     verts = []
    #     heights = []
    #     for row_idx, col_idx in zip(rows, cols):
    #         row = row_idx + y_start
    #         col = col_idx + x_start
    #         x = (col + 0.5) * res_x + offset_x
    #         y = bounds.top - (row + 0.5) * res_y + offset_y  # 修正y坐标
    #         z = 0
    #         verts.append((x, y, z))
    #         heights.append(float(data_chunk[row_idx, col_idx]))

    #     mesh = bpy.data.meshes.new(obj_name + "_mesh")
    #     mesh.from_pydata(verts, [], [])
    #     mesh.update()

    #     height_layer = mesh.attributes.new(name="height", type="FLOAT", domain="POINT")
    #     for i, height in enumerate(heights):
    #         height_layer.data[i].value = height

    #     point_cloud_obj = bpy.data.objects.new(obj_name + "_points", mesh)
    #     return point_cloud_obj

    def setup_geometry_nodes(self, obj, cell_size, props):
        # 创建新节点组并强制清理旧数据
        node_group_name = f"GN_RasterViz_{obj.name}"
        if node_group_name in bpy.data.node_groups:
            bpy.data.node_groups.remove(bpy.data.node_groups[node_group_name])

        node_group = bpy.data.node_groups.new(node_group_name, "GeometryNodeTree")
        interface = node_group.interface

        # 添加输入接口
        interface.new_socket(
            name="Geometry", in_out="INPUT", socket_type="NodeSocketGeometry"
        )
        height_scale = interface.new_socket(
            name="Height Scale", in_out="INPUT", socket_type="NodeSocketFloat"
        )
        height_scale.default_value = 1.0

        cube_size = interface.new_socket(
            name="Cube Size", in_out="INPUT", socket_type="NodeSocketFloat"
        )
        cube_size.default_value = cell_size

        # 添加输出接口
        interface.new_socket(
            name="Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry"
        )
        nodes = node_group.nodes
        links = node_group.links

        # 清空默认节点
        for node in nodes:
            nodes.remove(node)

        # 创建节点链 --------------------------------------------------------------
        ## 获取高度属性
        attribute = nodes.new("GeometryNodeInputNamedAttribute")
        attribute.data_type = "FLOAT"
        attribute.inputs["Name"].default_value = "height"
        attribute.location = (-800, -200)

        ## 数学运算
        # --- ABSOLUTE VALUE ---
        absolute_value = nodes.new(type="ShaderNodeMath")
        absolute_value.location = (-400, -500)
        absolute_value.operation = "ABSOLUTE"
        # --- MULTIPLY (1) ---
        multiply_1 = nodes.new(type="ShaderNodeMath")
        multiply_1.location = (-600, -300)
        multiply_1.operation = "MULTIPLY"
        multiply_1.inputs[1].default_value = 0.5
        # --- SYMBOL ---
        symbol = nodes.new(type="ShaderNodeMath")
        symbol.location = (-400, -300)
        symbol.operation = "SIGN"
        # --- MULTIPLY (2) ---
        multiply_2 = nodes.new(type="ShaderNodeMath")
        multiply_2.location = (-200, -300)
        multiply_2.operation = "MULTIPLY"
        multiply_2.inputs[1].default_value = 0.5

        # translation_math = nodes.new("ShaderNodeMath")
        # translation_math.operation = "MULTIPLY"
        # translation_math.label = "translation_z"
        # translation_math.inputs[1].default_value = 0.5  # 中心调整

        # scale_math = nodes.new("ShaderNodeMath")
        # scale_math.operation = "DIVIDE"
        # scale_math.label = "scale_z"

        # temp_multiply = nodes.new("ShaderNodeMath")
        # temp_multiply.operation = "MULTIPLY"

        # translation_vector = nodes.new("ShaderNodeCombineXYZ")
        # scale_vector = nodes.new("ShaderNodeCombineXYZ")
        # scale_vector.inputs["X"].default_value = 1.0
        # scale_vector.inputs["Y"].default_value = 1.0

        # --- COMBINE XYZ (1) ---
        combine_xyz_1 = nodes.new(type="ShaderNodeCombineXYZ")
        combine_xyz_1.location = (0, -300)
        combine_xyz_1.inputs[2].default_value = 1.0

        # 实例化系统
        cube = nodes.new("GeometryNodeMeshCube")
        cube.location = (-800, 200)

        # --- COMBINE XYZ (2) ---
        combine_xyz_2 = nodes.new(type="ShaderNodeCombineXYZ")
        combine_xyz_2.location = (-600, 200)
        combine_xyz_2.inputs[0].default_value = 1.0
        combine_xyz_2.inputs[1].default_value = 1.0

        # --- COMBINE XYZ (3) ---
        combine_xyz_3 = nodes.new(type="ShaderNodeCombineXYZ")
        combine_xyz_3.location = (-400, 200)
        combine_xyz_3.inputs[0].default_value = 0
        combine_xyz_3.inputs[1].default_value = 0

        instance = nodes.new("GeometryNodeInstanceOnPoints")
        instance.location = (-600, 0)

        # --- TRANSLATE INSTANCE ---
        translate_instance = nodes.new(type="GeometryNodeTranslateInstances")
        translate_instance.location = (-100, 0)
        # --- SCALE INSTANCE ---
        scale_instance = nodes.new(type="GeometryNodeScaleInstances")
        scale_instance.location = (150, 0)
        scale_instance.inputs["Scale"].default_value = (1, 1, 1)

        realize = nodes.new("GeometryNodeRealizeInstances")
        realize.location = (-200, 0)

        group_input = nodes.new("NodeGroupInput")
        group_input.location = (-1200, 0)

        group_output = nodes.new("NodeGroupOutput")
        group_output.location = (600, 0)

        # 正确连接节点 ------------------------------------------------------------
        # node_group.links.new(
        #     group_input.outputs["Geometry"], attribute.inputs[0]
        # )  # Input Geometry
        links = node_group.links
        # ---geometry
        links.new(group_input.outputs["Geometry"], instance.inputs["Points"])
        links.new(cube.outputs["Mesh"], instance.inputs["Instance"])
        links.new(instance.outputs["Instances"], translate_instance.inputs["Instances"])
        links.new(
            translate_instance.outputs["Instances"], scale_instance.inputs["Instances"]
        )
        links.new(scale_instance.outputs["Instances"], realize.inputs["Geometry"])
        links.new(realize.outputs["Geometry"], group_output.inputs["Geometry"])

        # ---height & translate & scale
        links.new(attribute.outputs["Attribute"], multiply_1.inputs[0])
        links.new(group_input.outputs["Height Scale"], multiply_1.inputs[1])
        links.new(multiply_1.outputs["Value"], symbol.inputs[0])
        links.new(multiply_1.outputs["Value"], absolute_value.inputs[0])
        links.new(symbol.outputs["Value"], multiply_2.inputs[0])
        links.new(multiply_2.outputs["Value"], combine_xyz_3.inputs["Z"])
        links.new(
            combine_xyz_3.outputs["Vector"], translate_instance.inputs["Translation"]
        )
        links.new(absolute_value.outputs["Value"], combine_xyz_2.inputs["Z"])
        links.new(combine_xyz_2.outputs["Vector"], scale_instance.inputs["Scale"])

        links.new(group_input.outputs["Cube Size"], combine_xyz_1.inputs["X"])
        links.new(group_input.outputs["Cube Size"], combine_xyz_1.inputs["Y"])
        links.new(combine_xyz_1.outputs["Vector"], cube.inputs["Size"])

        # node_group.links.new(attribute.outputs["Attribute"], translation_math.inputs[0])
        # node_group.links.new(
        #     group_input.outputs["Height Scale"], translation_math.inputs[1]
        # )  # Height Scale
        # node_group.links.new(
        #     translation_math.outputs[0], translation_vector.inputs["Z"]
        # )

        # node_group.links.new(attribute.outputs["Attribute"], temp_multiply.inputs[0])
        # node_group.links.new(
        #     group_input.outputs["Height Scale"], temp_multiply.inputs[1]
        # )  # Height Scale
        # node_group.links.new(temp_multiply.outputs[0], scale_math.inputs[0])
        # node_group.links.new(
        #     group_input.outputs["Cube Size"], scale_math.inputs[1]
        # )  # Cube Size
        # node_group.links.new(scale_math.outputs[0], scale_vector.inputs["Z"])

        # node_group.links.new(
        #     group_input.outputs["Cube Size"], cube.inputs["Size"]
        # )  # Cube Size

        # node_group.links.new(group_input.outputs["Geometry"], instance.inputs["Points"])
        # node_group.links.new(cube.outputs["Mesh"], instance.inputs["Instance"])

        # node_group.links.new(
        #     instance.outputs["Instances"], transform.inputs["Geometry"]
        # )
        # node_group.links.new(
        #     translation_vector.outputs["Vector"], transform.inputs["Translation"]
        # )
        # node_group.links.new(scale_vector.outputs["Vector"], transform.inputs["Scale"])

        # node_group.links.new(transform.outputs["Geometry"], realize.inputs["Geometry"])
        # node_group.links.new(
        #     realize.outputs["Geometry"], group_output.inputs["Geometry"]
        # )  # Output Geometry
        # # 位置调整（可选）
        # attribute.location = (-600, -200)
        # translation_math.location = (-400, -200)
        # translation_vector.location = (-200, -200)
        # scale_math.location = (-400, -400)
        # temp_multiply.location = (-600, -400)
        # scale_vector.location = (-200, -400)
        # cube.location = (-400, 200)
        # instance.location = (-200, 0)
        # transform.location = (0, 0)
        # realize.location = (200, 0)
        # 应用修改器 --------------------------------------------------------------
        gn_mod = obj.modifiers.new("GeometryNodes", "NODES")
        gn_mod.node_group = node_group
        gn_mod["Input_2"] = obj  # 设置输入几何体为点云数据
        gn_mod["Input_3"] = props.height_scale
        gn_mod["Input_4"] = cell_size


#    def setup_geometry_nodes(self, obj, data, nodata, props, res_x, res_y, width, height):
#        """为对象设置几何节点栅格系统"""
#        # 创建几何节点组
#        node_group_name = f"GN_RasterViz_{obj.name}"
#        if node_group_name in bpy.data.node_groups:
#            node_group = bpy.data.node_groups[node_group_name]
#        else:
#            node_group = bpy.data.node_groups.new(node_group_name, 'GeometryNodeTree')
#
#            # 设置输入输出
#            group_in = node_group.nodes.new('NodeGroupInput')
#            group_in.location = (-800, 0)
#            group_out = node_group.nodes.new('NodeGroupOutput')
#            group_out.location = (600, 0)
#
#            # 添加输入输出接口
#            node_group.inputs.new('NodeSocketFloat', 'Height Scale')
#            node_group.inputs.new('NodeSocketFloat', 'Cube Size')
#            node_group.outputs.new('NodeSocketGeometry', 'Geometry')
#
#            # 创建网格网格节点
#            grid = node_group.nodes.new('GeometryNodeMeshGrid')
#            grid.location = (-600, 0)
#            grid.inputs['Size X'].default_value = width * res_x
#            grid.inputs['Size Y'].default_value = height * res_y
#            grid.inputs['Vertices X'].default_value = width
#            grid.inputs['Vertices Y'].default_value = height
#
#            # 创建网格数据传输节点
#            attribute = node_group.nodes.new('GeometryNodeInputNamedAttribute')
#            attribute.location = (-600, -200)
#            attribute.data_type = 'FLOAT'
#            attribute.inputs['Name'].default_value = "Height"
#
#            # 创建立方体实例
#            cube = node_group.nodes.new('GeometryNodeMeshCube')
#            cube.location = (-400, 200)
#
#            # 创建实例化节点
#            instance = node_group.nodes.new('GeometryNodeInstanceOnPoints')
#            instance.location = (-200, 0)
#
#            # 创建高度转换节点
#            height_math = node_group.nodes.new('ShaderNodeMath')
#            height_math.location = (-400, -200)
#            height_math.operation = 'MULTIPLY'
#
#            # 创建变换节点
#            transform = node_group.nodes.new('GeometryNodeTransform')
#            transform.location = (0, 0)
#
#            # 创建实例化节点
#            realize = node_group.nodes.new('GeometryNodeRealizeInstances')
#            realize.location = (200, 0)
#
#            # 创建删除组件节点 (用于基于高度值筛选)
#            delete = node_group.nodes.new('GeometryNodeDeleteGeometry')
#            delete.location = (400, 0)
#            delete.domain = 'INSTANCE'
#
#            # 连接节点
#            node_group.links.new(grid.outputs['Mesh'], instance.inputs['Points'])
#            node_group.links.new(cube.outputs['Mesh'], instance.inputs['Instance'])
#            node_group.links.new(group_in.outputs['Height Scale'], height_math.inputs[1])
#            node_group.links.new(attribute.outputs['Attribute'], height_math.inputs[0])
#            node_group.links.new(height_math.outputs[0], transform.inputs['Translation'].driver_add('z').driver)
#            node_group.links.new(group_in.outputs['Cube Size'], cube.inputs['Size'])
#            node_group.links.new(instance.outputs['Instances'], transform.inputs['Geometry'])
#            node_group.links.new(transform.outputs['Geometry'], realize.inputs['Geometry'])
#            node_group.links.new(realize.outputs['Geometry'], group_out.inputs['Geometry'])
#
#        # 在对象上应用几何节点修改器
#        gn_mod = obj.modifiers.new("GeometryNodes", 'NODES')
#        gn_mod.node_group = node_group
#
#        # 设置节点组输入值
#        gn_mod["Input_1"] = props.raster_scale  # Height Scale
#        gn_mod["Input_2"] = props.cube_size * (1 - props.gap_size)  # Cube Size
#
#        # 创建高度属性数据
#        # 将栅格数据转换为高度属性
#        heights = data.flatten()
#        valid_mask = heights != nodata
#        heights = heights.astype(np.float32)
#
#        # 创建高度属性
#        height_attr_name = f"height_data_{obj.name}"
#        if height_attr_name in bpy.data.attributes:
#            height_attr = bpy.data.attributes[height_attr_name]
#        else:
#            height_attr = bpy.data.attributes.new(height_attr_name, 'FLOAT', 'POINT')
#            height_attr.data.values = heights[valid_mask]
#
#        # 将高度数据应用到几何节点
#        # 注意：此处需要进一步实现，在Blender API中可能需要使用驱动器或自定义脚本


class GIS_OT_LoadAttributes(Operator):
    bl_idname = "gis.load_attributes"
    bl_label = "加载属性字段"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        props = context.scene.gis_props
        if not props.shp_path or not os.path.exists(props.shp_path):
            self.report({"ERROR"}, "请选择有效的SHP文件")
            return {"CANCELLED"}

        try:
            # 读取SHP文件并获取属性字段
            gdf = gpd.read_file(props.shp_path)
            field_names = [col for col in gdf.columns if col != "geometry"]

            # 创建临时操作符来显示字段选择对话框
            # 由于Blender API限制，这里简化为直接打印字段名称
            self.report({"INFO"}, f"可用字段: {', '.join(field_names)}")

            return {"FINISHED"}
        except Exception as e:
            self.report({"ERROR"}, f"加载属性字段失败: {str(e)}")
            return {"CANCELLED"}


class GIS_OT_OptimizeMemory(Operator):
    bl_idname = "gis.optimize_memory"
    bl_label = "优化内存使用"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        try:
            # 收集垃圾
            bpy.ops.outliner.orphans_purge(do_recursive=True)
            self.report({"INFO"}, "内存优化完成")
            return {"FINISHED"}
        except Exception as e:
            self.report({"ERROR"}, f"内存优化失败: {str(e)}")
            return {"CANCELLED"}


class GIS_PT_Panel(Panel):
    bl_label = "GIS 可视化"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "工具"

    def draw(self, context):
        layout = self.layout
        props = context.scene.gis_props

        # 矢量数据部分
        box = layout.box()
        box.label(text="矢量数据")
        row = box.row()
        row.prop(props, "shp_path")

        row = box.row()
        row.prop(props, "name_field")
        row.prop(props, "height_field")

        row = box.row()
        row.prop(props, "height_scale")
        row.prop(props, "height_offset")

        row = box.row()
        row.operator("gis.load_attributes", text="加载属性字段")
        row.operator("gis.import_vector", text="导入矢量")

        # 栅格数据部分
        box = layout.box()
        box.label(text="栅格数据")
        row = box.row()
        row.prop(props, "tif_path")

        row = box.row()
        row.prop(props, "raster_scale")

        # row = box.row()
        # row.prop(props, "cube_size")
        row.prop(props, "gap_size")

        # row = box.row()
        # row.prop(props, "chunk_size")

        row = box.row()
        row.operator("gis.import_raster", text="导入栅格")

        # 高级选项
        box = layout.box()
        box.label(text="高级选项")

        row = box.row()
        row.prop(props, "auto_center")

        row = box.row()
        row.prop(props, "use_lod")
        row.prop(props, "lod_distance")

        # row = box.row()
        # row.prop(props, "color_map")

        # row = box.row()
        # row.prop(props, "render_mode")

        row = box.row()
        row.operator("gis.optimize_memory", text="优化内存")


def register():
    bpy.utils.register_class(GisProps)
    bpy.utils.register_class(GIS_OT_ImportVector)
    bpy.utils.register_class(GIS_OT_ImportRaster)
    bpy.utils.register_class(GIS_OT_LoadAttributes)
    bpy.utils.register_class(GIS_OT_OptimizeMemory)
    bpy.utils.register_class(GIS_PT_Panel)
    bpy.types.Scene.gis_props = PointerProperty(type=GisProps)


def unregister():
    bpy.utils.unregister_class(GIS_PT_Panel)
    bpy.utils.unregister_class(GIS_OT_OptimizeMemory)
    bpy.utils.unregister_class(GIS_OT_LoadAttributes)
    bpy.utils.unregister_class(GIS_OT_ImportRaster)
    bpy.utils.unregister_class(GIS_OT_ImportVector)
    bpy.utils.unregister_class(GisProps)
    del bpy.types.Scene.gis_props


if __name__ == "__main__":
    register()
