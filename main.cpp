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

#include "DualGraph.h"

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
	GraphNode *graphNode;
	Point center;
	std::string colorStr;
	std::set<Vertex*> vertices;
	std::set<Face*> faces;
	int pid;
	std::vector<Halfedge*> boundary;
	std::vector<Vertex*> corners;
	std::vector<Vertex*> sampledPoints;
};


void generateSeeds(Mesh *mesh, std::vector<Patch*> &patches, int N) {
	char buf[512];

	for (int i = 0; i < N; ++i) {
		
		//int vid = rand() % mesh->numVertices() + 1;
		int vid = mesh->numVertices() % N + i;

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
	int pid = 0;
	for (auto p : patches) {
		std::cout << "tracing patch: " << ++pid << "\n";

		traceBoundary(p);
		//output boundary
		std::set<Vertex *> vertices;
		for (auto he : p->boundary) {
			vertices.insert(he->source());
		}
	}
}


double avg_length = 0;
DualGraph* generateGraph(std::vector<Patch*>& patches) {
	
	DualGraph *dualGraph = new DualGraph;

	std::map<Vertex*, std::vector<Patch*> > ver2patchMap;
	for (auto p : patches) {

		GraphNode *node = dualGraph->addNode();
		p->graphNode = node;

		for (auto v : p->vertices) {
			auto pvec = ver2patchMap[v];
			for (size_t i = 0; i < pvec.size(); ++i) {
				Patch *p0 = pvec[i];
				dualGraph->addEdge(p0->graphNode, node);
				dualGraph->addEdge(node, p0->graphNode);
			}
			ver2patchMap[v].push_back(p);
		}
	}


	//well.. let's print out the dualgraph

	//convert the dual graph to cm
	//first, the position, each center of the cluster
	std::ofstream dgfile("dualgraph.cm");
	
	for (auto p : patches) {
		GraphNode *gnode = p->graphNode;
		dgfile << "Vertex " << gnode->id() + 1 << " " << p->center[0] << " " << p->center[1] << " " << p->center[2] << "\n";
	}

	for (auto p : patches) {
		GraphNode *gnode = p->graphNode;
		for (GraphNode::EdgeIter eit(gnode); !eit.end(); ++eit) {
			GraphEdge *edge = *eit;
			if (edge->to()->id() > gnode->id()) {
				dgfile << "Edge " << edge->from()->id() + 1 << " " << edge->to()->id() + 1 << "\n";
			}
		}
	}

	dualGraph->saveMetis("dualgraph.metis");

	dgfile.close();
	


	std::vector <std::string> edgestr;
	std::set<Vertex*> vertices;
	double total_boundary_length = 0;
	int boundary_count = 0;
	for (auto p : patches) {
		std::set<Vertex*> pvertices;
		for (auto he : p->boundary) {
			
			if (!he->twin()) {
				total_boundary_length += (he->source()->point() - he->target()->point()).norm();
				++boundary_count;
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

	avg_length = total_boundary_length / boundary_count;
	std::map<Vertex*, int> prevOrder, nowOrder;
	std::ofstream output("graph.m");

	int vid = 1;
	for (auto v : vertices) {
		
		output << "Vertex " << vid << " " << v->point()[0] << " " << v->point()[1] << " " << v->point()[2] << "\n";
		prevOrder[v] = v->index();
		nowOrder[v] = vid++;
		
	}

	int fid = 0;
	for (auto p : patches) {
		output << "Face " << ++fid << " ";
	
		for (auto v : p->corners) {
			
			output << nowOrder[v]  << " ";
		}
		output << "\n";
	}

	for (auto v : vertices) {
		
		v->index() = prevOrder[v];
	}

	return dualGraph;
}

void sampling (std::vector<Patch*>& patches, double threhold) {
	//collect all the edges..
	//typedef std::pair<Vertex*, Vertex*> VertexPair;
	//std::set< VertexPair > pairset;
	typedef std::pair< std::pair<Vertex*, Vertex*>, int> SampledPoint;
	std::map<SampledPoint, Vertex*> spmap;
	for (auto p : patches) {
		std::cout << "patch!\n";
		size_t cornerSize = p->corners.size();
		for (size_t i = 0 ; i < cornerSize; ++i) {
			Vertex *v0 = p->corners[i];
			Vertex *v1 = p->corners[(i + 1) % cornerSize];
		
			std::cout <<"v0: " << v0->point()[0] << " " << v0->point()[1] << " " << v0->point()[2] << "\n";
			std::cout <<"v1: " << v1->point()[0] << " " << v1->point()[1] << " " << v1->point()[2] << "\n";

			const Point p0 = v0->point();
			const Point p1 = v1->point();
			double length = (p0 - p1).norm();

			p->sampledPoints.push_back(v0);

			if (length > threhold) {
				int numpt = length / threhold;
		
				for (int j = 1; j < numpt; ++j) {
					SampledPoint sp = std::make_pair( std::make_pair(v0, v1), j);
					//well, looking for the vertex is it's been created.
					Vertex *v_new = spmap[sp];
					if (!v_new) {
						double lambda = (double) j / numpt;
						Point newp = p0 * (1- lambda) + p1 * lambda;
						v_new = new Vertex;
						v_new->index() = -1;
						v_new->point() = newp;

						spmap[sp] = v_new;
						sp = std::make_pair(std::make_pair(v1, v0), numpt - j);
						spmap[sp] = v_new;
					} else {
						std::cout << "found match!  ";
					}

					std::cout << "  new pt" <<  v_new->point()[0] << " " <<v_new->point()[1] << " " << v_new->point()[2] << "\n";
					for (int j = 0; j < p->sampledPoints.size(); ++j) {
						Point chkp = p->sampledPoints[j]->point();
						if ( (chkp - v_new->point()).norm() < 1e-6) {
							std::cout << "cat!\n";
						}
					}
					p->sampledPoints.push_back(v_new);
				}
			} 

		
		}
	}

	//okay, let's save out all the patches..

	char buf[128];
	
	for (size_t i = 0; i < patches.size(); ++i) {
		Patch *p = patches[i];

		sprintf_s(buf, "patch_%02d", i + 1);
		std::ofstream output(buf);
		int vid = 0;
		size_t spSize = p->sampledPoints.size();
		for (size_t j = 0; j < spSize; ++j) {
			Vertex *v = p->sampledPoints[j];
			
			output << "Vertex " << ++vid << " " << v->point()[0] << " " << v->point()[1] << " " << v->point()[2];
			if (v->index() < 0)  {
				output << " {new}";
			}
			output << "\n";
		}

		
		for (size_t j = 0; j < spSize; ++j) {
			output << "Edge " <<  j + 1 << " " << (j + 1) % spSize + 1 << "\n";
		}

		output.close();
	}


	std::ofstream sampled("sampled.m");
	std::set<Vertex*> verset;
	std::map<Vertex*, int> nowOrder;
	int vid = 0;
	for (auto p : patches) {
		for (auto v: p->sampledPoints) {
			auto it = verset.insert(v);
			if (it.second) {

				nowOrder[v] = ++vid;
				sampled << "Vertex " << vid << " " << v->point()[0] << " " << v->point()[1] << " " << v->point()[2] << "\n";
			}
		}
	}


	int fid = 0;
	for (auto p : patches) {
		sampled << "Face " << ++fid << " ";
		for (auto v : p->sampledPoints) {
			sampled << nowOrder[v] << " ";
		}
		sampled << "\n";
	}



	sampled.close();


}

int main(int argc, char *argv[]) {
	if (argc != 4) { 
		std::cout << argv[0] << "input.m patchnum sampling\n";
		exit(-1);
	}
	
	int patchnum = atoi(argv[2]);
	double threshold = atof(argv[3]);
	srand (static_cast <unsigned> (time(0)));

	//load seeds;
	std::vector<Patch*> patches;
	//loadSeeds("seed_sim.m", patches);
	//loadSeeds("seed2.m", patches);

	Mesh *mesh = new Mesh;
	mesh->readMFile(argv[1]);

	generateSeeds(mesh, patches, patchnum);
	
	for (int i = 0; i <= 500; ++i) {
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
	DualGraph *dualGraph = generateGraph(patches);

	sampling(patches, avg_length);

	delete dualGraph;
	delete mesh;
	return 0;
} 