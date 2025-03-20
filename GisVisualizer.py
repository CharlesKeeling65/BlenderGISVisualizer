# -*- coding: utf-8 -*-
"""
@File    :   GisVisualizer.py
@Time    :   2023/06/15
@Author  :   Charles Keeling (Wang Yubo)
@Version :   1.0
@Desc    :   用于Blender中导入并可视化GIS数据的插件，支持矢量数据和栅格数据的处理和展示
"""
import math
import os
import tempfile
import time
from functools import lru_cache
from typing import Any, Dict, List, Optional, Tuple, Union

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

# 插件信息
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
    """GIS数据属性组，存储用户界面配置信息"""
    
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
    """导入矢量数据的操作器类"""
    bl_idname = "gis.import_vector"
    bl_label = "导入矢量数据"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context: bpy.types.Context) -> Dict[str, str]:
        """执行导入矢量数据的操作
        
        Args:
            context: Blender的上下文对象
            
        Returns:
            Dict[str, str]: 操作状态字典
        """
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

    def get_field_names(self, gdf: gpd.GeoDataFrame) -> List[str]:
        """获取并返回GeoDataFrame的所有字段名称
        
        Args:
            gdf: GeoDataFrame数据对象
            
        Returns:
            List[str]: 字段名称列表
        """
        return list(gdf.columns)

    def process_vector(self, context: bpy.types.Context, props: GisProps) -> None:
        """处理矢量数据，导入到Blender场景中
        
        Args:
            context: Blender的上下文对象
            props: GIS属性配置
        """
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

    def create_polygon(
        self, 
        collection: bpy.types.Collection, 
        geom: Union[Polygon, MultiPolygon], 
        row: gpd.GeoSeries, 
        idx: Union[int, str], 
        props: GisProps, 
        offset: Vector
    ) -> None:
        """根据几何数据创建Blender多边形对象
        
        Args:
            collection: 目标集合对象
            geom: 几何数据对象
            row: GeoDataFrame的一行数据
            idx: 对象索引
            props: GIS属性配置
            offset: 坐标偏移量
        """
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
    """导入栅格数据的操作器类"""
    bl_idname = "gis.import_raster"
    bl_label = "导入栅格数据"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context: bpy.types.Context) -> Dict[str, str]:
        """执行导入栅格数据的操作
        
        Args:
            context: Blender的上下文对象
            
        Returns:
            Dict[str, str]: 操作状态字典
        """
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

    def process_raster(self, context: bpy.types.Context, props: GisProps) -> None:
        """处理栅格数据，导入到Blender场景中
        
        Args:
            context: Blender的上下文对象
            props: GIS属性配置
        """
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


    def setup_geometry_nodes(
        self, 
        obj: bpy.types.Object, 
        cell_size: float, 
        props: GisProps
    ) -> None:
        """设置几何节点系统用于可视化栅格数据
        
        Args:
            obj: 栅格点云对象
            cell_size: 栅格单元大小
            props: GIS属性配置
        """
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
        translate_instance.inputs['Local Space'].default_value = False
        # --- SCALE INSTANCE ---
        scale_instance = nodes.new(type="GeometryNodeScaleInstances")
        scale_instance.location = (150, 0)
        scale_instance.inputs["Scale"].default_value = (1, 1, 1)
        # 关闭缩放实例的局部空间选项
        scale_instance.inputs['Local Space'].default_value = False

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

        # 应用修改器 --------------------------------------------------------------
        gn_mod = obj.modifiers.new("GeometryNodes", "NODES")
        gn_mod.node_group = node_group
        gn_mod["Input_2"] = obj  # 设置输入几何体为点云数据
        gn_mod["Input_3"] = props.raster_scale
        gn_mod["Input_4"] = cell_size


class GIS_OT_LoadAttributes(Operator):
    """加载属性字段的操作器类"""
    bl_idname = "gis.load_attributes"
    bl_label = "加载属性字段"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context: bpy.types.Context) -> Dict[str, str]:
        """执行加载属性字段操作
        
        Args:
            context: Blender的上下文对象
            
        Returns:
            Dict[str, str]: 操作状态字典
        """
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
    """优化内存使用的操作器类"""
    bl_idname = "gis.optimize_memory"
    bl_label = "优化内存使用"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context: bpy.types.Context) -> Dict[str, str]:
        """执行内存优化操作
        
        Args:
            context: Blender的上下文对象
            
        Returns:
            Dict[str, str]: 操作状态字典
        """
        try:
            # 收集垃圾
            bpy.ops.outliner.orphans_purge(do_recursive=True)
            self.report({"INFO"}, "内存优化完成")
            return {"FINISHED"}
        except Exception as e:
            self.report({"ERROR"}, f"内存优化失败: {str(e)}")
            return {"CANCELLED"}


class GIS_PT_Panel(Panel):
    """GIS可视化工具面板类"""
    bl_label = "GIS 可视化"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "工具"

    def draw(self, context: bpy.types.Context) -> None:
        """绘制面板UI
        
        Args:
            context: Blender的上下文对象
        """
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


def register() -> None:
    """注册插件类和属性"""
    bpy.utils.register_class(GisProps)
    bpy.utils.register_class(GIS_OT_ImportVector)
    bpy.utils.register_class(GIS_OT_ImportRaster)
    bpy.utils.register_class(GIS_OT_LoadAttributes)
    bpy.utils.register_class(GIS_OT_OptimizeMemory)
    bpy.utils.register_class(GIS_PT_Panel)
    bpy.types.Scene.gis_props = PointerProperty(type=GisProps)


def unregister() -> None:
    """注销插件类和属性"""
    bpy.utils.unregister_class(GIS_PT_Panel)
    bpy.utils.unregister_class(GIS_OT_OptimizeMemory)
    bpy.utils.unregister_class(GIS_OT_LoadAttributes)
    bpy.utils.unregister_class(GIS_OT_ImportRaster)
    bpy.utils.unregister_class(GIS_OT_ImportVector)
    bpy.utils.unregister_class(GisProps)
    del bpy.types.Scene.gis_props


if __name__ == "__main__":
    register()
