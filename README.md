# BlenderGISVisualizer (Blender GIS 可视化插件)

This is a Blender plugin for 3D visualization of vector or raster map files.

这是一个 Blender 插件，用于将矢量或栅格类型的地图文件进行 3D 可视化

---

## 效果示意图

栅格图
![栅格图](./asset/raster_map.png)
矢量图
![矢量图](./asset/vector_map.png)

---

## 使用方法

### 1. 在 Blender 中安装依赖

在 Blender 对应的 Python 环境中安装依赖

```bash
<blender_python_interpreter> -m pip install -r requirements.txt
```

### 2. 在 Blender 中安装插件

- 在 Blender 中打开 `Edit` -> `Preferences` -> `Add-ons` -> `Install...`
- 选择本项目中的 `GisVisualizer.py` 文件
- 点击 `Install Add-on`
- 在插件列表中找到 `Blender GIS Visualizer`，勾选启用

---

## 任务列表

### Done

- [x] 实现了.shp .gpkg .geojson 等矢量格式地图的导入和绘制
- [x] 实现了.tif 栅格格式地图的导入和绘制
- [x] 矢量格式地图可以进行统一以高度调整渐变颜色

### Todo

- [ ] 准备继续用几何节点实现栅格 3d 图的颜色渐变材质
- [ ] 准备支持矢量图和栅格图地理坐标系的变换
