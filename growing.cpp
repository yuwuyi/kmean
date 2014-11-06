#include "growing.h"


void RegionGrowing::initPatches(Face *face) {
	GrowingPatch *p = new GrowingPatch;
	p->center = calcCenter(face);
	for (FaceHalfedgeIterator fheit(face); !fheit.end(); ++fheit) {
		p->frontEdges.push_back(*fheit);
	}
}

bool RegionGrowing::grow() {
	for (size_t i = 0; i < m_patches.size(); ++i) {
		GrowingPatch *patch = m_patches[i];
	}

	return true;
}

bool RegionGrowing::grow(GrowingPatch *patch) {

	std::set<Face*, GrowingPatch> candidateFaces(*patch);
	std::set<Halfedge*> candidateFrontHalfEdges;
	for (size_t i = 0; i < patch->frontEdges.size(); ++i) {
		Halfedge *he = patch->frontEdges[i];
		Face *f0 = he->face();
		Halfedge *ohe = he->twin();
		if (!ohe) {
			continue;
		}
		Face *f1 = ohe->face();

		if (face2patchMap[f0]) {
			candidateFrontHalfEdges.insert(he);
			candidateFaces.insert(f0);
		} 
	
		if (face2patchMap[f1]) {
			candidateFrontHalfEdges.insert(ohe);
			candidateFaces.insert(f1);
		} 
	}

	//update the map
	for (auto f : candidateFaces) {
		face2patchMap[f] = patch;
		f->PropertyStr() =  "rgb=(1 0 0)";
	}

	if (!candidateFrontHalfEdges.empty()) {
		patch->frontEdges.clear();
		//update the front edge;
		for (auto he : candidateFrontHalfEdges) {
			Face *f0 = he->face();
			Halfedge *ohe = he->twin();
			if (ohe) {
				Face *f1 = ohe->face();
				if (!face2patchMap[f0]) {
					patch->frontEdges.push_back(he);
				}
				if (!face2patchMap[f1]) {
					patch->frontEdges.push_back(ohe);
				}
			}
		}
	}

	return true;
}