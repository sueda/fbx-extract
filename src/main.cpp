#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <set>
#define _USE_MATH_DEFINES
#include <cmath>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>

#include "ofbx.h"

using namespace std;
using namespace glm;

// https://github.com/nem0/OpenFBX/blob/master/demo/main.cpp
// https://github.com/Larry955/FbxParser/tree/master/FbxParser

string FILENAME;

class MyVertex
{
public:
	MyVertex()
	{
		p.x = p.y = p.z = t.x = t.y = n.x = n.y = n.z = 0.0f;
	}
	
	bool operator==(const MyVertex &v) const
	{
		float thresh = 1e-6;
		auto vec3Eq = [thresh](const vec3 &a, const vec3 &b)
		{
			return fabs(a.x - b.x) < thresh && fabs(a.y - b.y) < thresh && fabs(a.z - b.z) < thresh;
		};
		auto vec2Eq = [thresh](const vec2 &a, const vec2 &b)
		{
			return fabs(a.x - b.x) < thresh && fabs(a.y - b.y) < thresh;
		};
		if(!vec3Eq(p, v.p)) {
			return false;
		}
		if(!vec2Eq(t, v.t)) {
			return false;
		}
		if(!vec3Eq(n, v.n)) {
			return false;
		}
		return true;
	}
	
	vec3 p;
	vec2 t;
	vec3 n;
	vector<float> w;
	vector<int> i;
};

class MyTriangle
{
public:
	MyTriangle()
	{
		v[0] = v[1] = v[2] = 0;
	}
	int v[3];
};

class MyMesh
{
public:
	string name;
	vector<MyVertex> verts;
	vector<MyTriangle> tris;
	
	// We want to store the indices of unique vertices, since there
	// are many duplicated vertices.
	// Let's say that verts contains duplicate vertices A and D:
	//                      0 1 2 3 4 5 6 7
	//    verts          = [A B a C D a d E] <- duplicates in lower case
	//    vertsUnique    = [0 1 3 4 7]
	//    vertsUniqueMap = [0 1 0 2 3 0 3 4] <- indexes into vertsUnique
	vector<int> vertsUnique;
	vector<int> vertsUniqueMap;
	
	int maxInfluences; // maximum number of bone influences
};

vector<MyMesh> myMeshes;

bool saveGeom(const ofbx::IScene *scene)
{
	cout << "=== geometry ===" << endl;
	int mesh_count = scene->getMeshCount();
	myMeshes.resize(mesh_count);
	for (int k = 0; k < mesh_count; ++k) {
		const ofbx::Mesh *mesh = scene->getMesh(k);
		const ofbx::Geometry *geom = mesh->getGeometry();
		int vertex_count = geom->getVertexCount();
		int index_count = geom->getIndexCount();
		cout << mesh->name << ": " << vertex_count << " verts, ";
		myMeshes[k].name = mesh->name;
		replace(myMeshes[k].name.begin(), myMeshes[k].name.end(), ':', '_');
		replace(myMeshes[k].name.begin(), myMeshes[k].name.end(), '/', '_');
		
		myMeshes[k].verts.resize(vertex_count);
		auto &verts = myMeshes[k].verts;
		const ofbx::Vec3* vertices = geom->getVertices();
		for (int i = 0; i < vertex_count; ++i) {
			const ofbx::Vec3 &v = vertices[i];
			verts[i].p.x = v.x;
			verts[i].p.y = v.y;
			verts[i].p.z = v.z;
		}

		bool has_normals = geom->getNormals() != nullptr;
		if(has_normals) {
			const ofbx::Vec3* normals = geom->getNormals();
			// This will fail if ofbx::LoadFlags::TRIANGULATE is not used
			assert(geom->getIndexCount() == vertex_count);
			for(int i = 0; i < vertex_count; ++i) {
				const ofbx::Vec3 &n = normals[i];
				verts[i].n.x = n.x;
				verts[i].n.y = n.y;
				verts[i].n.z = n.z;
			}
		}

		bool has_uvs = geom->getUVs() != nullptr;
		if(has_uvs) {
			const ofbx::Vec2 *uvs = geom->getUVs();
			// This will fail if ofbx::LoadFlags::TRIANGULATE is not used
			assert(geom->getIndexCount() == vertex_count);
			for (int i = 0; i < vertex_count; ++i) {
				const ofbx::Vec2 &uv = uvs[i];
				verts[i].t.x = uv.x;
				verts[i].t.y = uv.y;
			}
		}

		// Assume that we are working with a triangulated array, so
		// there is no index/element array.
		assert(vertex_count == index_count);
		int tri_count = vertex_count/3;
		myMeshes[k].tris.resize(tri_count);
		auto &tris = myMeshes[k].tris;
		for(int i = 0; i < tri_count; ++i) {
			// This is index-0. Convert to index-1 before exporting to OBJ.
			tris[i].v[0] = 3*i + 0;
			tris[i].v[1] = 3*i + 1;
			tris[i].v[2] = 3*i + 2;
		}

		// Find duplicates (slow O(n^2))
		// E.g.,
		//                      0 1 2 3 4 5 6 7
		//    verts          = [A B a C D a d E] <- duplicates in lower case
		//    vertsUnique    = [0 1 3 4 7]
		//    vertsUniqueMap = [0 1 0 2 3 0 3 4] <- indexes into vertsUnique
		auto &vertsUnique = myMeshes[k].vertsUnique;
		auto &vertsUniqueMap = myMeshes[k].vertsUniqueMap;
		vertsUniqueMap.resize(vertex_count, -1);
		for(int i = 0; i < vertex_count; ++i) {
			int duplj = -1;
			for(int j = 0; j < (int)vertsUnique.size(); ++j) {
				if(verts[i] == verts[vertsUnique[j]]) {
					// Found duplicate in vertsUnique
					duplj = j;
					break;
				}
			}
			if(duplj == -1) {
				// No duplicate was found, so create a new entry
				vertsUniqueMap[i] = vertsUnique.size();
				vertsUnique.push_back(i);
			} else {
				// Duplicate found
				vertsUniqueMap[i] = duplj;
			}
		}
		cout << vertsUnique.size() << " unique" << endl;
	}
	
	// Save
	for(auto &mesh : myMeshes) {
		string path = FILENAME + "_" + mesh.name + ".obj";
		FILE* fp = fopen(path.c_str(), "wb");
		if (!fp) {
			cout << "Could not write to " << path << endl;
			return false;
		}
		cout << "Saving to " << path << endl;
		fprintf(fp, "# %s\n", mesh.name.c_str());
		const auto &verts = mesh.verts;
		const auto &tris = mesh.tris;
		const auto &vertsUnique = mesh.vertsUnique;
		const auto &vertsUniqueMap = mesh.vertsUniqueMap;
		fprintf(fp, "# %lu vertices, %lu triangles\n", vertsUnique.size(), tris.size());
		for(int i : vertsUnique) {
			const MyVertex &v = verts[i];
			fprintf(fp, "v %f %f %f\n", v.p.x, v.p.y, v.p.z);
		}
		for(int i : vertsUnique) {
			const MyVertex &v = verts[i];
			fprintf(fp, "vt %f %f\n", v.t.x, v.t.y);
		}
		for(int i : vertsUnique) {
			const MyVertex &v = verts[i];
			fprintf(fp, "vn %f %f %f\n", v.n.x, v.n.y, v.n.z);
		}
		for(const MyTriangle &t : mesh.tris) {
			// Use vertsUniqueMap to index into the vector of unique verts.
			// Add 1 to convert to index-1 from index-0.
			int v0 = vertsUniqueMap[t.v[0]] + 1;
			int v1 = vertsUniqueMap[t.v[1]] + 1;
			int v2 = vertsUniqueMap[t.v[2]] + 1;
			fprintf(fp, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", v0, v0, v0, v1, v1, v1, v2, v2, v2);
		}
		fclose(fp);
	}
	
	return true;
}

bool saveSkin(const ofbx::IScene *scene)
{
	cout << "=== skin ===" << endl;
	
	// Get sizes
	int mesh_count = scene->getMeshCount();
	assert(mesh_count == myMeshes.size());
	int cluster_count = 0;
	for(int i = 0; i < mesh_count; ++i) {
		const ofbx::Mesh *mesh = scene->getMesh(i);
		const ofbx::Geometry *geom = mesh->getGeometry();
		const ofbx::Skin *skin = geom->getSkin();
		cout << mesh->name << ": ";
		int vertex_count = geom->getVertexCount();
		cout << vertex_count << " vertices" << endl;
		int cluster_count_ = skin->getClusterCount();

		// All meshes have the same number of clusters
		if(i == 0) {
			cluster_count = cluster_count_;
		} else {
			assert(cluster_count_ == cluster_count);
		}
	}

	// Populate attachment info
	for(int i = 0; i < mesh_count; ++i) {
		const ofbx::Mesh *mesh = scene->getMesh(i);
		const ofbx::Geometry *geom = mesh->getGeometry();
		const ofbx::Skin *skin = geom->getSkin();
		myMeshes[i].maxInfluences = 0;
		auto &verts = myMeshes[i].verts;
		for(int j = 0; j < cluster_count; ++j) {
			const ofbx::Cluster *cluster = skin->getCluster(j);
			int indices_count = cluster->getIndicesCount();
			int weights_count = cluster->getWeightsCount();
			assert(indices_count == weights_count);
			if(weights_count > 0) {
				const int *indices = cluster->getIndices();
				const double *weights = cluster->getWeights();
				for(int k = 0; k < weights_count; ++k) {
					int index = indices[k];
					double weight = weights[k];
					assert(index < geom->getVertexCount());
					assert(0.0 <= weight && weight <= 1.0);
					verts[index].w.push_back(weight);
					verts[index].i.push_back(j);
					if(verts[index].i.size() > myMeshes[i].maxInfluences) {
						myMeshes[i].maxInfluences = (int)verts[index].i.size();
					}
				}
			}
		}
	}
	
	// Write to attachment file
	for(const auto &mesh : myMeshes) {
		string path = FILENAME + "_" + mesh.name + "_skin.txt";
		FILE* fp = fopen(path.c_str(), "wb");
		if (!fp) {
			cout << "Could not write to " << path << endl;
			return false;
		}
		const auto &verts = mesh.verts;
		const auto &vertsUnique = mesh.vertsUnique;
		int vertex_count = vertsUnique.size();
		cout << "Saving to " << path << endl;
		fprintf(fp, "# 1st line: vertCount boneCount maxInfluences\n");
		fprintf(fp, "# Each subsequent line corresponds to a vertex.\n");
		fprintf(fp, "# In each line, the first number is the number of bone influences for the vertex.\n");
		fprintf(fp, "# The next numbers are \"influence\" pairs of {bone index, bone weight}.\n");
		fprintf(fp, "%d %d %d\n", vertex_count, cluster_count, mesh.maxInfluences);
		for(int i : vertsUnique) {
			const MyVertex &v = verts[i];
			int influences = v.w.size();
			assert(influences == v.i.size());
			fprintf(fp, "%d ", influences);
			for(int j = 0; j < influences; ++j) {
				fprintf(fp, "%d %f ", v.i[j], v.w[j]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
	
	return true;
}

void traverseAll(const ofbx::Object *object, int depth)
{
	int depthMax = 6; // -1 to disable
	if(depth == depthMax) {
		return;
	}
	string s;
	for(int i = 0; i < depth; ++i) {
		s += "\t";
	}
	const char* typeLabel;
	switch (object->getType()) {
		case ofbx::Object::Type::GEOMETRY: typeLabel = "GEOMETRY"; break;
		case ofbx::Object::Type::MESH: typeLabel = "MESH"; break;
		case ofbx::Object::Type::MATERIAL: typeLabel = "MATERIAL"; break;
		case ofbx::Object::Type::ROOT: typeLabel = "ROOT"; break;
		case ofbx::Object::Type::TEXTURE: typeLabel = "TEXTURE"; break;
		case ofbx::Object::Type::NULL_NODE: typeLabel = "NULL"; break;
		case ofbx::Object::Type::LIMB_NODE: typeLabel = "LIMB"; break;
		case ofbx::Object::Type::NODE_ATTRIBUTE: typeLabel = "ATTRIBUTE"; break;
		case ofbx::Object::Type::CLUSTER: typeLabel = "CLUSTER"; break;
		case ofbx::Object::Type::SKIN: typeLabel = "SKIN"; break;
		case ofbx::Object::Type::ANIMATION_STACK: typeLabel = "ANIM_STACK"; break;
		case ofbx::Object::Type::ANIMATION_LAYER: typeLabel = "ANIM_LAYER"; break;
		case ofbx::Object::Type::ANIMATION_CURVE: typeLabel = "ANIM_CURVE"; break;
		case ofbx::Object::Type::ANIMATION_CURVE_NODE: typeLabel = "ANIM_CURVE_NODE"; break;
		default: assert(false); break;
	}
	
	if(object->getType() != ofbx::Object::Type::NODE_ATTRIBUTE) {
		cout << s << typeLabel << " : " << object->name;
		if(object->getType() == ofbx::Object::Type::GEOMETRY) {
			//cout << " : " << ((ofbx::Geometry*)object)->getVertexCount() << " verts";
		}
		cout << endl;
	}
	
	int i = 0;
	while(const ofbx::Object* child = object->resolveObjectLink(i)) {
		traverseAll(child, depth + 1);
		++i;
	}
}

void traverseLimbs(const ofbx::Object *limb, vector<const ofbx::Object*> &limbVec, map<const ofbx::Object*, int> &limbMap, vector<int> &limbParents)
{
	if(limb->getType() != ofbx::Object::Type::LIMB_NODE) {
		return;
	}
	
	limbVec.push_back(limb);
	limbMap[limb] = limbMap.size();
	const ofbx::Object *parent = limb->getParent();
	if(parent != nullptr && parent->getType() == ofbx::Object::Type::LIMB_NODE) {
		limbParents.push_back(limbMap[parent]);
	} else {
		limbParents.push_back(-1);
	}
	
	int i = 0;
	while(const ofbx::Object* child = limb->resolveObjectLink(i)) {
		traverseLimbs(child, limbVec, limbMap, limbParents);
		++i;
	}
}

bool saveAnim(const ofbx::IScene *scene)
{
	cout << "=== skeleton ===" << endl;
	
	// http://docs.autodesk.com/FBX/2014/ENU/FBX-SDK-Documentation/index.html?url=cpp_ref/class_fbx_pose.html,topicNumber=cpp_ref_class_fbx_pose_html3b4fd8fc-688e-43ec-b1df-56acb1cce550
	// Cluster is a bone
	// Use depthMax = 6 in traverseAll() to see the hierarchy:
	//
	//	MESH : newVegas:Elvis_BodyGeo
	//		GEOMETRY : newVegas:Elvis_BodyGeo
	//			SKIN : Skin newVegas:Elvis_BodyGeo
	//				CLUSTER : Link 0
	//					LIMB : newVegas:Hips
	//					LIMB : newVegas:Hips
	//				CLUSTER : Link 1
	//					LIMB : newVegas:Pelvis
	//					LIMB : newVegas:Pelvis
	//				CLUSTER : Link 2
	//					LIMB : newVegas:LeftUpLeg
	//					LIMB : newVegas:LeftUpLeg
	//				CLUSTER : Link 3
	//					LIMB : newVegas:LeftLeg
	//					LIMB : newVegas:LeftLeg
	//
	// TransformLink is the global transform of the bone(link) at the binding moment
	
	// Assume that each mesh skin has the same number of clusters (bones) as all the
	// other mesh skins.
	int mesh_count = scene->getMeshCount();
	if(mesh_count == 0) {
		cout << "This file has no mesh" << endl;
		return true;
	}
	
	// Casts double[16] to glm::mat4
	auto toMat4 = [](const double *array) {
		mat4 M;
		for(int row = 0; row < 4; ++row) {
			for(int col = 0; col < 4; ++col) {
				double v = array[4 * col + row];
				M[col][row] = v;
			}
		}
		return M;
	};
	
	// Get an Euler angle rotation matrix
	auto toR = [](const ofbx::Vec3 &v, ofbx::RotationOrder ro) {
		mat4 I = mat4(1.0);
		float x = v.x * M_PI / 180.0f;
		float y = v.y * M_PI / 180.0f;
		float z = v.z * M_PI / 180.0f;
		mat4 Rx = rotate(I, x, vec3(1.0, 0.0, 0.0));
		mat4 Ry = rotate(I, y, vec3(0.0, 1.0, 0.0));
		mat4 Rz = rotate(I, z, vec3(0.0, 0.0, 1.0));
		mat4 R = I;
		switch(ro) {
			case ofbx::RotationOrder::EULER_XYZ:
				R = Rz * Ry * Rx;
				break;
			case ofbx::RotationOrder::EULER_XZY:
				R = Ry * Rz * Rx;
				break;
			case ofbx::RotationOrder::EULER_YZX:
				R = Rx * Rz * Ry;
				break;
			case ofbx::RotationOrder::EULER_YXZ:
				R = Rz * Rx * Ry;
				break;
			case ofbx::RotationOrder::EULER_ZXY:
				R = Ry * Rx * Rz;
				break;
			case ofbx::RotationOrder::EULER_ZYX:
				R = Rx * Ry * Rz;
				break;
			case ofbx::RotationOrder::SPHERIC_XYZ:
				assert(false);
				break;
		}
		return R;
	};
	
	// Get a translation matrix
	auto toT = [](const ofbx::Vec3 &v) {
		mat4 I = mat4(1.0);
		mat4 T = translate(I, vec3(v.x, v.y, v.z));
		return T;
	};
	
	// Get a scale matrix
	auto toS = [](const ofbx::Vec3 &v) {
		mat4 I = mat4(1.0);
		if(v.x == 0 && v.y == 0 && v.z == 0) {
			// Assuming that 0 scale was sent in by mistake
			return I;
		} else {
			mat4 S = scale(I, vec3(v.x, v.y, v.z));
			return S;
		}
	};
	
	// prints glm::mat4
	auto printMat4 = [](const mat4 &M, const char *name = nullptr) {
		if(name != nullptr) {
			cout << name << " = [" << endl;
		}
		for(int row = 0; row < 4; ++row) {
			for(int col = 0; col < 4; ++col) {
				cout << M[col][row] << " ";
			}
			cout << endl;
		}
		if(name != nullptr) {
			cout << "];" << endl;
		}
	};
	
	// Find the root limb (hips)
	const ofbx::Object *root = nullptr; // root limb
	int i = 0;
	while(const ofbx::Object* child = scene->getRoot()->resolveObjectLink(i)) {
		if(child->getType() == ofbx::Object::Type::LIMB_NODE) {
			root = child;
			break;
		}
		++i;
	}
	
	// Traverse hierarchy
	vector<const ofbx::Object*> limbVec;
	map<const ofbx::Object*, int> limbMap;
	vector<int> limbParents;
	traverseLimbs(root, limbVec, limbMap, limbParents);
	assert(limbVec.size() == limbParents.size());
	
	// Find the mesh with the maximum number of clusters
	int cluster_count = 0;
	const ofbx::Skin *skin = nullptr;
	for(int i = 0; i < mesh_count; ++i) {
		const ofbx::Mesh *mesh = scene->getMesh(i);
		const ofbx::Geometry *geom = mesh->getGeometry();
		if(geom->getSkin()->getClusterCount() > cluster_count) {
			skin = geom->getSkin();
			cluster_count = skin->getClusterCount();
			break;
		}
	}
	assert(skin != nullptr);
	cout << cluster_count << " clusters" << endl;
	
	// Store rest pose
	int key_count_max = -1;
	vector<mat4> pose0;
	pose0.reserve(cluster_count);
	for(int j = 0; j < cluster_count; ++j) {
		const ofbx::Cluster *cluster = skin->getCluster(j);
		//cout << "cluster " << j << ", " << cluster->name << endl;
		const ofbx::Matrix m = cluster->getTransformLinkMatrix();
		// Need to cast from double to float
		// Assume both are column major
		pose0.push_back(toMat4(m.m));
		// Get key_count, used later
		const ofbx::Object *limb = cluster->getLink();
		if(limb->getType() == ofbx::Object::Type::LIMB_NODE) {
			if(root == nullptr) {
				root = limb;
			}
			// limb has one or two children that are curve nodes.
			int i = 0;
			while(const ofbx::Object* child = limb->resolveObjectLink(i)) {
				if(child->getType() == ofbx::Object::Type::ANIMATION_CURVE_NODE) {
					const ofbx::AnimationCurveNode* node = (ofbx::AnimationCurveNode*)child;
					const ofbx::AnimationCurve* curveX = node->getCurve(0);
					if(curveX->getKeyCount() > key_count_max) {
						key_count_max = curveX->getKeyCount();
					}
				}
				++i;
			}
		}
	}
	assert(key_count_max != -1);
	
	// Store the local transform values
	vector< vector<mat4> > Ts(cluster_count);
	vector< vector<mat4> > Rs(cluster_count);
	for(int j = 0; j < cluster_count; ++j) {
		const ofbx::Cluster *cluster = skin->getCluster(j);
		const ofbx::Object *limb = cluster->getLink();
		// limb has one or more children that are curve nodes.
		vector<ofbx::Vec3> rotations;
		vector<ofbx::Vec3> positions;
		int i = 0;
		int key_count = 0; // key count for this limb
		while(const ofbx::Object* child = limb->resolveObjectLink(i)) {
			if(child->getType() == ofbx::Object::Type::ANIMATION_CURVE_NODE) {
				const ofbx::AnimationCurveNode* node = (ofbx::AnimationCurveNode*)child;
				const ofbx::AnimationCurve* curveX = node->getCurve(0);
				const ofbx::AnimationCurve* curveY = node->getCurve(1);
				const ofbx::AnimationCurve* curveZ = node->getCurve(2);
				key_count = curveX->getKeyCount();
				assert(key_count == curveY->getKeyCount());
				assert(key_count == curveZ->getKeyCount());
				const float *xvals = curveX->getKeyValue();
				const float *yvals = curveY->getKeyValue();
				const float *zvals = curveZ->getKeyValue();
				for(int k = 0; k < key_count; ++k) {
					ofbx::Vec3 v;
					v.x = xvals[k];
					v.y = yvals[k];
					v.z = zvals[k];
					if(strcmp(child->name, "R") == 0) {
						rotations.push_back(v);
					} else if(strcmp(child->name, "T") == 0) {
						positions.push_back(v);
					} else {
						cout << "Unsupported animation curve type" << endl;
						assert(false);
					}
				}
			}
			++i;
		}
		//cout << limb->name << ": " << key_count << " frames" << endl;
		for(int k = 0; k < key_count; ++k) {
			// http://docs.autodesk.com/FBX/2014/ENU/FBX-SDK-Documentation/cpp_ref/class_fbx_cluster.html
			if(k < positions.size()) {
				Ts[j].push_back(toT(positions[k]));
			}
			if(k < rotations.size()) {
				Rs[j].push_back(toR(rotations[k], limb->getRotationOrder()));
			}
		}
	}
	
	// These poses are in local coordinates, so transform them to world
	// https://help.autodesk.com/view/FBX/2017/ENU/?guid=__files_GUID_10CDD63C_79C1_4F2D_BB28_AD2BE65A02ED_htm
	// WorldTransform = ParentWorldTransform * T * Roff * Rp * Rpre * R * Rpost^-1 * Rp^-1 * Soff * Sp * S * Sp^-1
	vector< vector<mat4> > Es(cluster_count);
	for(int k = 0; k < key_count_max; ++k) {
		for(int j = 0; j < cluster_count; ++j) {
			const ofbx::Object *limb = limbVec[j];
			int p = limbParents[j];
			mat4 P = mat4(1.0);
			if(p != -1) {
				P = Es[p][k];
			}
			ofbx::RotationOrder ro = limb->getRotationOrder();
			mat4 R = toR(limb->getLocalRotation(), ro);
			mat4 T = toT(limb->getLocalTranslation());
			if(Rs[j].empty()) {
				// No key frames for this limb. Use rest pose.
			} else {
				if(k < Rs[j].size()) {
					// Use this key frame.
					R = Rs[j][k];
				} else {
					// Use the last stored key frame.
					R = Rs[j][Rs[j].size() - 1];
				}
			}
			if(Ts[j].empty()) {
				// No key frames for this limb. Use rest pose.
			} else {
				if(k < Ts[j].size()) {
					// Use this key frame.
					T = Ts[j][k];
				} else {
					// Use the last stored key frame.
					T = Ts[j][Ts[j].size() - 1];
				}
			}
			mat4 Roff = toR(limb->getRotationOffset(), ro);
			mat4 Rp = toR(limb->getRotationPivot(), ro);
			mat4 Rpre = toR(limb->getPreRotation(), ro);
			mat4 Rpost = toR(limb->getPostRotation(), ro);
			mat4 Soff = toS(limb->getScalingOffset());
			mat4 Sp = toS(limb->getScalingPivot());
			mat4 S = toS(limb->getLocalScaling());
			//printMat4(P, "P");
			//printMat4(R, "R");
			//printMat4(T, "T");
			//printMat4(Roff, "Roff");
			//printMat4(Rp, "Rp");
			//printMat4(Rpost, "Rpost");
			//printMat4(Soff, "Soff");
			//printMat4(Sp, "Sp");
			//printMat4(S, "S");
			mat4 E = P * T * Roff * Rp * Rpre * R * inverse(Rpost) * inverse(Rp) * Soff * Sp * S * inverse(Sp);
			Es[j].push_back(E);
		}
	}
	
	// Write to skeleton file
	string path = FILENAME + "_skel.txt";
	FILE* fp = fopen(path.c_str(), "wb");
	if (!fp) {
		cout << "Could not write to " << path << endl;
		return false;
	}
	cout << "Saving to " << path << endl;
	fprintf(fp, "# 1st line: frameCount boneCount\n");
	fprintf(fp, "# Each subsequent line is a frame with \"boneCount\" sets of\n");
	fprintf(fp, "# quaternions (x,y,z,w) and positions (x,y,z).\n");
	fprintf(fp, "%d %d\n", key_count_max, cluster_count);
	
	for(const auto &M : pose0) {
		quat q = quat_cast(M);
		vec3 p(M[3]);
		fprintf(fp, "%f %f %f %f %f %f %f ", q.x, q.y, q.z, q.w, p.x, p.y, p.z);
	}
	fprintf(fp, "\n");
	for(int k = 0; k < key_count_max; ++k) {
		for(int j = 0; j < cluster_count; ++j) {
			const mat4 &E = Es[j][k];
			quat q = quat_cast(E);
			vec3 p(E[3]);
			fprintf(fp, "%f %f %f %f %f %f %f ", q.x, q.y, q.z, q.w, p.x, p.y, p.z);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
	return true;
}

// https://www.oreilly.com/library/view/c-cookbook/0596007612/ch10s15.html
string getFileName(const string& s)
{
	char sep = '/';
#ifdef _WIN32
	sep = '\\';
#endif
	size_t i = s.rfind(sep, s.length());
	if(i != string::npos) {
		return(s.substr(i+1, s.length() - i));
	}
	return("");
}

int main(int argc, char **argv)
{
	if(argc < 2) {
		cout << "Usage: fbx-extract <FBX_FILE>" << endl;
		return -1;
	}
	
	// Get file handle
	string filename = argv[1];
	FILE* fp = fopen(filename.c_str(), "rb");
	if(!fp) {
		cout << filename << " not found" << endl;
		return -1;
	}

	// Get file size and allocate memory
	fseek(fp, 0, SEEK_END);
	long file_size = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	auto* content = new ofbx::u8[file_size];
	
	// Load data into memory
	fread(content, 1, file_size, fp);
	
	// Parse
	cout << "Parsing " << filename << endl;
	const ofbx::IScene* scene = ofbx::load((ofbx::u8*)content, file_size, (ofbx::u64)ofbx::LoadFlags::TRIANGULATE);
	if(!scene) {
		cout << ofbx::getError() << endl;
	}
	
	// Extract just the filename
	FILENAME = getFileName(filename);
	FILENAME = FILENAME.substr(0, FILENAME.length() - 4);
	
	// DEBUG: Traverse tree
	//traverseAll(scene->getRoot(), 0);
	
	// Parse and save
	saveGeom(scene);
	saveSkin(scene);
	saveAnim(scene);
	
	// Delete data
	delete [] content;
	
	cout << "done" << endl;
}
