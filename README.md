# fbx-extract

This program parses an FBX file and writes six types of output files.

1. An OBJ file of the mesh
2. An ASCII file with skinning information
3. An ASCII file with the skeletal animation information
4. Texture files
5. An ASCII file with the input information
6. Four ASCII files with the local transforms for each joint in skeletal animation information

The input FBX file can be downloaded from [Mixamo](mixamo.com).

Authors:
- Shinjiro Sueda
- Ziyan Xiong


## Building and running

This project builds on top of the following projects:

- [GLM](https://glm.g-truc.net/0.9.9/)
- [OpenFBX](https://github.com/nem0/OpenFBX) (included)
- [miniz](https://github.com/richgel999/miniz) (included)

To build the project, run `cmake` and `make`, and then run with the FBX file as an argument:

```
> mkdir build
> cd build
> cmake ..
> make
> ./fbx-extract ../data/bigvegas_Walking.fbx
```


## Output files

### Geometry

- A standard OBJ file with vertex positions, normals, and texture coordinates.
- Faces are triangulated.
- Common vertices are merged.


### Skinning weights

- Filename ends with `_skin.txt`.
- Comments start with `#`.
- The first line has three integers:
  - Vertex count: should be the same as the corresponding OBJ file
  - Bone count: the number of bones in the animation
  - Max influences: this is the maximum number of non-zero weights per vertex
- Each subsequent line corresponds to a vertex.
  - In each line, the first number is the number of bone influences for the vertex.
  - The next numbers are "influence" pairs of {bone index, bone weight}. The bone index is 0-indexed.

```
4583 82 9
1 20 1.000000 
1 20 1.000000 
2 19 0.037664 20 0.962336 
```


### Skeletal animation

- Filename ends with `_skel.txt`.
- Comments start with `#`.
- The first line has two integers:
  - Frame count: the number of frames in the animation
  - Bone count: the number of bones in the animation
- Each subsequent line is a frame with "boneCount" sets of quaternions `(x,y,z,w)` and positions `(x,y,z)`.
- The first frame is for the rest pose, so there are actually `frameCount + 1` lines.

```
27 82
0.000000 0.000000 0.000000 1.000000 0.000000 96.330055 9.825774 ...
0.008630 0.039924 0.024404 0.998867 0.657645 87.535782 2.391427 ...
```


### Textures

- All of the textures contained in the input FBX file, usually are PNG or JPG files.


### Input file

- The input information for this [Assignment](https://people.engr.tamu.edu/sueda/courses/CSCE489/2021F/assignments/A2/index.html). 
- This file only lists the name of the first extracted texture. If it is not the ideal one, please edit it manually.


### Local transforms

- There are four ASCII files storing the local transforms information, so we can convert them back without loss.
- The world transform can be restored using: `WorldTransform = ParentWorldTransform * T * Roff * Rp * Rpre * R * Rpost^-1 * Rp^-1 * Soff * Sp * S * Sp^-1`, described [here](https://help.autodesk.com/view/FBX/2017/ENU/?guid=__files_GUID_10CDD63C_79C1_4F2D_BB28_AD2BE65A02ED_htm).
- Comments start with '#'.

1. `*_static_transforms.txt` contains the static transforms between the joints.
  - The first line contains the number of joints.
  - Each line next contains 8 matrices: `T Roff Rp Rpre Rpost Soff Sp S`.
  - Each matrix has 7 numbers: 4 for quaternion `(x, y, z, w)` and 3 for position `(x, y, z)`, so that each line has `7 * 8 = 56` numbers.
  
```
83
0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 ...
0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 ...
0.000000 0.000000 0.000000 1.000000 18.735231 -5.666687 -4.541196 ...
```

2. `*_skel_local.txt` contains the Euler angles of local transformation for each joint in each frame.
  - The first line has two integers: `frameCount` and `boneCount`, which are the same as the skeletal animation file.
  - Each subsequent line is a frame with `boneCount` sets of Euler angles.
  - The rotation order of the Euler angles is stored in the next file.
  - The root joint contains Euler angles and positions, while all other joints contain just Euler angles. There will be `3 + 3 * boneCount numbers` on each line.
  
```
27 83
0.019251 0.079420 0.049618 0.657645 87.535782 2.391427 ...
0.016953 0.083084 0.038091 0.214333 87.434326 2.976235 ...
```

3. `*_hierarchy.txt` contains the parent-child hierarchy for each joint in each frame.
  - The first line contains the number of joints.
  - Each subsequent line is in the form of `<JOINT INDEX> <PARENT INDEX> <ROTATION ORDER> <JOINT NAME>`
  - If a joint has `-1` as its parent index, then it is the root joint.

```
83
0 -1 EULER_XYZ newVegas:Hips
1 0 EULER_XYZ newVegas:Pelvis
```

4. `_binding_pose_local.txt` contains the Euler angle of local transformation for each joint of the binding pose.
  - It has the same format as the `_skel_local.txt` file.

```
27 83
0.000000 0.000000 0.000000 0.000000 96.330055 9.825774...
```


## TODO

- Add support for using FBX files with only skeletal animation data. Currently, the FBX file must contain geometry, skinning, and animation data. This means that to download new animations for an existing character, you need to download another FBX file with redundant information.
