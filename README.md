# Offsets and Remeshing

This repository contains an implementation of the method presented in "Feature-Preserving Offset Mesh Generation from Topology-Adapted Octrees".
It generates offsets for triangle meshes and performs feature aware remeshing.
Both components, offsetting, and remeshing can be used independently.

```
offr::offr
```

## Dependencies

The code is based on [CGAL](https://www.cgal.org/) and was tested with version 5.5.2. For the offsetting library requires slight modification of the Octrees package.

> TODO list required modifications

## Offsetting

This library implements topology aware octree refinement and Dual Contouring.

```
offr::dc
```

## Remeshing

This library implements the remeshing. Although it was specifically designed for offset remeshing, we extended the implementation to work with any type of hermite data.

```
offr::remesh
```
