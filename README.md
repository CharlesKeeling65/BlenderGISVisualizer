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
