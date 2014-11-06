#include <set>
#include <map>
#include <random>
#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "MeshLib_Core/Point.h"
#include "MeshLib_Core/Mesh.h"
#include "MeshLib_Core/Iterators.h"

#include "growing.h"

Point center(Face *face) {
	Point point;
	for (FaceVertexIterator fvit(face); !fvit.end(); ++fvit) {
		Vertex *v = *fvit;
		point += v->point();
	}

	point /= 3.0;
	return point;
}

double distance_linfinity(const Point & p1, const Point& p2) {
	Point p = p1 - p2;
	p[0] = fabs(p[0]);
	p[1] = fabs(p[1]);
	p[2] = fabs(p[2]);
	return *std::max_element(p.v, p.v + 2);
}

struct Patch {
	Point center;
	std::string colorStr;
	std::set<Vertex*> vertices;
	std::set<Face*> faces;
	int pid;
	std::vector<Halfedge*> boundary;
	std::vector<Vertex*> corners;
};


void generateSeeds(Mesh *mesh, std::vector<Patch*> &patches, int N) {
	char buf[512];

	for (int i = 0; i < N; ++i) {
		
		int vid = rand() % mesh->numVertices() + 1;

		Patch *patch = new Patch;
		patch->center = mesh->indVertex(vid)->point();

		float r = static_cast <float> (rand() % 100) / static_cast <float> (100);
		float g = static_cast <float> (rand() % 100) / static_cast <float> (100);
		float b = static_cast <float> (rand() % 100) / static_cast <float> (100);


		sprintf_s(buf, "rgb=(%f %f %f)", r, g, b);
		patch->colorStr = buf;
		
		patches.push_back(patch);
	}
}

void loadSeeds(const char *filename, std::vector<Patch*> &patches) {
	std::ifstream input(filename);
	std::string line;
	char buf[512];
	while(input.good()) {
		getline(input, line);
		if(line.empty()) {
			continue;
		}
		std::stringstream ss(line);
		Point pt;
		ss >> pt.v[0] >> pt.v[1] >> pt.v[2];
		Patch *patch = new Patch;
		patch->center = pt;

		float r = static_cast <float> (rand() % 100) / static_cast <float> (100);
		float g = static_cast <float> (rand()% 100) / static_cast <float> (100);
		float b = static_cast <float> (rand()% 100) / static_cast <float> (100);
		 
		sprintf_s(buf, "rgb=(%f %f %f)", r, g, b);
		patch->colorStr = buf;
			
		

		patches.push_back(patch);
	}
	input.close();
}

std::set<Face*> clustered_faces;
Patch* closestPatch(const Point& point, const std::vector<Patch*>& patches) {
	double mindist = 1e10;
	Patch *resultPatch = nullptr;
	for (auto &p : patches) {
		double dist = distance_linfinity(point, p->center);
		if (dist < mindist) {
			mindist = dist;
			resultPatch = p;
		}
	}
	return resultPatch;
}


void checkPatches(const std::vector<Patch*>& patches) {
	int i = 1;
	for (auto &p : patches) {
		
		if (p->vertices.empty()) {
			std::cout << i << " get empty faces!\n";
		}
		++i;
	}
}

void clustering(Mesh *mesh,  const	std::vector<Patch*> & patches) {

	for (auto p : patches) {
		p->vertices.clear();
		p->faces.clear();
	}

	for (MeshFaceIterator mfit(mesh); !mfit.end(); ++mfit) {
		Face *f = *mfit;
		//compute the center;
		Point c = center(f);
		
		
		Patch *p = closestPatch(c, patches);
	
		f->PropertyStr() = p->colorStr;
		for (FaceVertexIterator fvit(f); !fvit.end(); ++fvit) {
			p->vertices.insert(*fvit);
		}
		p->faces.insert(f);
	}
}

Point compute_center(Patch* patch) {
	Point p;
 	for (auto &v : patch->vertices) {
		p += v->point();
	}
	return p / patch->vertices.size();
}

void update(const std::vector<Patch*>& patches) {
	for (auto &p : patches) {
		p->center = compute_center(p);
		
	}
}

void saveCenters(const char *filename, const std::vector<Patch*>& patches) {
	std::ofstream output(filename);
	int i = 1;
	for (auto &p : patches) {
		output << "Vertex " << i << " " << p->center[0] << " " <<p->center[1] << " " << p->center[2] << "\n";
		++i;
	}
	output.close();
}
//
//int growing_main(int argc, char *argv[]) {
//	if (argc != 3) {
//		std::cout << argv[0] << "input.m output.m\n";
//		exit(-1);
//	}
//
//
//}

void traceBoundary(Patch *patch) {
	std::list<Halfedge*> hes;
	for (auto face : patch->faces) {
		for (FaceHalfedgeIterator fheit(face); !fheit.end(); ++fheit) {
			Halfedge *he = *fheit;
			Halfedge *ohe = he->twin();
			if(!ohe) {
				hes.push_back(he);
				continue;
			}

			Face *oface = ohe->face();
			if (patch->faces.find(oface) == patch->faces.end()) {
				hes.push_back(he);
			}
		}
	}

	Halfedge *start_e = hes.front();
	std::vector<Halfedge*>& boundary = patch->boundary;
	boundary.push_back(start_e);

	//boundary.erase(boundary.begin());
	while (!hes.empty()) {
		auto heit = hes.begin();
		bool isProgress = false;
		while (heit != hes.end())
		{
			Halfedge *next_e = *heit;
			if (start_e->target() == next_e->source()) { 
				boundary.push_back(next_e);
				hes.erase(heit++); 
				start_e = next_e;
				isProgress = true;
			}
			else
			{
				++heit;
			}
		}

		if (!isProgress) {
			std::cout << "\t\tfailed!\n";
			break;
		}
	}
	
}

void traceBoundary(std::vector<Patch*>& patches) {
	char buf[128];
	int pid = 0;
	for (auto p : patches) {
		std::cout << "tracing patch: " << ++pid << "\n";

		traceBoundary(p);
		//output boundary
		std::set<Vertex *> vertices;
		for (auto he : p->boundary) {
			vertices.insert(he->source());
		}
		
		sprintf_s(buf, "graph.%02d", pid);
		std::ofstream output(buf);
		for (auto v : vertices) {
			output << "Vertex " << v->index() + 1 << " " << v->point()[0] << " " << v->point()[1] << " " << v->point()[2] << "\n";
 		}
		for (auto he : p->boundary) {
			output << "Edge " << he->source()->index() + 1 << " " << he->target()->index() + 1 << "\n";
		}
		output.close();
		
		//sprintf_s(buf, "curve2m graph.%02d g%02d.m", pid, pid);
		//system(buf);
	}
}

void generateGraph(std::vector<Patch*>& patches) {
	std::map<Vertex*, std::vector<Patch*> > ver2patchMap;
	for (auto p : patches) {
		for (auto v : p->vertices) {
			ver2patchMap[v].push_back(p);
		}
	}

	std::vector <std::string> edgestr;
	std::set<Vertex*> vertices;
	for (auto p : patches) {
		std::set<Vertex*> pvertices;
		for (auto he : p->boundary) {
			
			if (!he->twin()) {
				if ( ver2patchMap[he->source()].size() == 1) {
					//continue;
				}
				vertices.insert(he->source());
				auto itpair = pvertices.insert(he->source());
				if (itpair.second) {
					p->corners.push_back(he->source());
				}

				if ( ver2patchMap[he->target()].size() == 1) {
				//	continue;
				}
				vertices.insert(he->target());
				itpair = pvertices.insert(he->target());
				if (itpair.second) {
					p->corners.push_back(he->target());
				}
				
			}  
			if ( ver2patchMap[he->source()].size() >= 3) {
				vertices.insert(he->source());
				auto itpair = pvertices.insert(he->source());			
				if (itpair.second) {
					p->corners.push_back(he->source());
				}
					
			}


			if ( ver2patchMap[he->target()].size() >= 3) {
				vertices.insert(he->target());
				auto itpair = pvertices.insert(he->target());			
				if (itpair.second) {
					p->corners.push_back(he->target());
				}

			}

		}
	}

	std::ofstream output("graph.cm");

	for (auto v : vertices) {
		output << "Vertex " << v->index() + 1<< " " << v->point()[0] << " " << v->point()[1] << " " << v->point()[2] << "\n";
	}


	/*for (auto p : patches) {
		size_t csize = p->corners.size();
		for (size_t i = 0; i < csize; ++i) {
			output << "Edge " << p->corners[i]->index() + 1 << " " << p->corners[(i + 1) % csize]->index() + 1 << "\n";
		}
	}
	*/

	int fid = 0;
	for (auto p : patches) {
		output << "Face " << ++fid << " ";
		
		
		for (auto v : p->corners) {
			
			output << v->index() + 1  << " ";
		}
		output << "\n";
	}
}

int main(int argc, char *argv[]) {
	if (argc != 3) { 
		std::cout << argv[0] << "input.m output.m\n";
		exit(-1);
	}
	
	srand (static_cast <unsigned> (time(0)));

	//load seeds;
	std::vector<Patch*> patches;
	//loadSeeds("seed_sim.m", patches);
	//loadSeeds("seed2.m", patches);

	Mesh *mesh = new Mesh;

	mesh->readMFile(argv[1]);

	generateSeeds(mesh, patches, 80);
	
	for (int i = 0; i <= 400; ++i) {
		/*if (i % 100 == 0) {
			sprintf_s(buf, "center_%d.cm", i+1);
			saveCenters(buf, patches);
		}*/
		clustering(mesh, patches);
		checkPatches(patches);
		update(patches);
		/*if (i % 100 == 0) {
			sprintf_s(buf, "out_%d.m", i+1);
			mesh->writeMFile(buf);
		}*/
		
	}
	mesh->writeMFile("end.m");
	std::cout << "end!\n";
	traceBoundary(patches);
	generateGraph(patches);
	delete mesh;
	return 0;
} 