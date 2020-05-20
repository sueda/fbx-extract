# fbx-extract
Extracts data from an FBX file

This program parses an FBX file and writes three types of output files.

1. An OBJ file of the mesh
2. An ascii file with skinning information
3. An ascii file with the skeletal animation information

The skinning file has the following format.
- Comments start with '#'.
- The first line has three integers:
  - Vertex count: should be the same as the corresponding OBJ file
  - Bone count: the number of bones in the animation
  - Max influences: this is the maximum number of non-zero weights per vertex
- Each subsequent line corresponds to a vertex.
  - In each line, the first number is the number of bone influences for the vertex.
  - The next numbers are "influence" pairs of {bone index, bone weight}.

```
4583 82 9
1 20 1.000000 
1 20 1.000000 
2 19 0.037664 20 0.962336 
```

The animation file has the following format.
- Comments start with '#'.
- The first line has two integers:
  - Frame count: the number of frames in the animation
  - Bone count: the number of bones in the animation
- Each subsequent line is a frame with "boneCount" sets of quaternions (x,y,z,w) and positions (x,y,z).
- The first frame is for the rest pose, so there are actually "frame count + 1" lines.

```
27 82
0.000000 0.000000 0.000000 1.000000 0.000000 96.330055 9.825774 ...
0.008630 0.039924 0.024404 0.998867 0.657645 87.535782 2.391427 ...
```

# Libraries used

- GLM <https://glm.g-truc.net/0.9.9/index.html>
- OpenFBX <https://github.com/nem0/OpenFBX>


# TODO

- Extract textures from the FBX file. Currently, you need to extract them manually using Blender.
