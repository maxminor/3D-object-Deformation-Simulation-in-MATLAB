//Computational Fabrication Assignment #1
#include <iostream>
#include <vector>
#include "../include/CompFab.h"
#include "../include/Mesh.h"

//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    /********* ASSIGNMENT *********/
    /* Ray-Triangle intersection test: Return 1 if ray intersects triangle, 
     * 0 otherwise */
     
     //Vector of edge of triangle
     CompFab::Vec3 ab = triangle.m_v2-triangle.m_v1;
     CompFab::Vec3 ac = triangle.m_v3-triangle.m_v1;
     CompFab::Vec3 bc = triangle.m_v3-triangle.m_v2;
     
     CompFab::Vec3 n = ab%ac;
     
     //in case the ray and the plane is parallel
     if (n*ray.m_direction == 0) return 0;
     
     //find coordinate of ray-plane intersection
     double d = triangle.m_v1*n;
     
     double t = (d-(n*ray.m_origin))/(n*ray.m_direction);
     //ray formula R(x) = P + t(d). t must not be negative
	 if (t < 0) return 0;
     
     CompFab::Vec3 dis;
     dis[0] = t*ray.m_direction[0];
     dis[1] = t*ray.m_direction[1];
     dis[2] = t*ray.m_direction[2];
     
     //check if ray-plane intersection point is inside of the triangle
     CompFab::Vec3 DoP = ray.m_origin+dis;
     
     CompFab::Vec3 ADoP = DoP-triangle.m_v1;
     CompFab::Vec3 BDoP = DoP-triangle.m_v2;
	 CompFab::Vec3 CDoP = DoP-triangle.m_v3;
	 
	 double fanA = (ab%ADoP)*n;
	 if (fanA < 0 ) return 0;
	 double fanB = (bc%BDoP)*n;
	 if (fanB < 0 ) return 0;
	 double fanC = (CDoP%ac)*n;
	 if (fanC < 0) return 0;

    
    return 1;


}

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxelPos, CompFab::Vec3 &dir)
{
    
    unsigned int numHits = 0;
    
    /********* ASSIGNMENT *********/
    /* Check and return the number of times a ray cast in direction dir, 
     * from voxel center voxelPos intersects the surface */
 	CompFab::Ray ray(voxelPos, dir);
 	for (int j = 0 ; j<g_triangleList.size(); j++) 
	 	numHits += (rayTriangleIntersection(ray, g_triangleList[j]));	
    
    return numHits;
}

bool loadMesh(char *filename, unsigned int dim)
{
    g_triangleList.clear();
    
    Mesh *tempMesh = new Mesh(filename, true);
    
    CompFab::Vec3 v1, v2, v3;

    //copy triangles to global list
    for(unsigned int tri =0; tri<tempMesh->t.size(); ++tri)
    {
        v1 = tempMesh->v[tempMesh->t[tri][0]];
        v2 = tempMesh->v[tempMesh->t[tri][1]];
        v3 = tempMesh->v[tempMesh->t[tri][2]];
        g_triangleList.push_back(CompFab::Triangle(v1,v2,v3));
    }

    //Create Voxel Grid
    CompFab::Vec3 bbMax, bbMin;
    BBox(*tempMesh, bbMin, bbMax);
    
    //Build Voxel Grid
    double bbX = bbMax[0] - bbMin[0];
    double bbY = bbMax[1] - bbMin[1];
    double bbZ = bbMax[2] - bbMin[2];
    double spacing;
    
    if(bbX > bbY && bbX > bbZ)
    {
        spacing = bbX/(double)(dim-2);
    } else if(bbY > bbX && bbY > bbZ) {
        spacing = bbY/(double)(dim-2);
    } else {
        spacing = bbZ/(double)(dim-2);
    }
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    g_voxelGrid = new CompFab::VoxelGrid(bbMin-hspacing, dim, dim, dim, spacing);

    delete tempMesh;
    
    return true;
   
}

void saveVoxelsToObj(const char * outfile)
{
 
    Mesh box;
    Mesh mout;
    int nx = g_voxelGrid->m_dimX;
    int ny = g_voxelGrid->m_dimY;
    int nz = g_voxelGrid->m_dimZ;
    double spacing = g_voxelGrid->m_spacing;
    
    CompFab::Vec3 hspacing(0.5*spacing, 0.5*spacing, 0.5*spacing);
    
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                if(!g_voxelGrid->isInside(ii,jj,kk)){
                    continue;
                }
                CompFab::Vec3 coord(0.5f + ((double)ii)*spacing, 0.5f + ((double)jj)*spacing, 0.5f+((double)kk)*spacing);
                CompFab::Vec3 box0 = coord - hspacing;
                CompFab::Vec3 box1 = coord + hspacing;
                makeCube(box, box0, box1);
                mout.append(box);
            }
        }
    }

    mout.save_obj(outfile);
}


int main(int argc, char **argv)
{

    unsigned int dim = 32; //dimension of voxel grid (e.g. 32x32x32)
    unsigned int big_dim = 64;
    unsigned int small_dim = 16;
    unsigned int even_small_dim = 8;

    //Load OBJ
    if(argc < 3)
    {
        std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename \n";
        return 0;
    }
    
    std::cout<<"Load Mesh : "<<argv[1]<<"\n";
    loadMesh(argv[1], big_dim);
    std::cout<<"done loading mesh" << std::endl;
    

    
    //Cast ray, check if voxel is inside or outside
    //even number of surface intersections = outside (OUT then IN then OUT)
    // odd number = inside (IN then OUT)
    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction(1.0,0.0,0.0);
    
    /********* ASSIGNMENT *********/
    /* Iterate over all voxels in g_voxelGrid and test whether they are inside our outside of the
     * surface defined by the triangles in g_triangleList */
    for (unsigned int i=0; i<g_voxelGrid->m_dimX; i++) {
    	for (unsigned int j= 0; j<g_voxelGrid->m_dimY; j++) {
    		for (unsigned int k= 0; k<g_voxelGrid->m_dimZ; k++) {
    			CompFab::Vec3 coor(g_voxelGrid->m_lowerLeft[0] + ((double)i)*g_voxelGrid->m_spacing,
								   g_voxelGrid->m_lowerLeft[1] + ((double)j)*g_voxelGrid->m_spacing, 
								   g_voxelGrid->m_lowerLeft[2] + ((double)k)*g_voxelGrid->m_spacing);
				g_voxelGrid->isInside(i, j, k) = (numSurfaceIntersections(coor, direction) % 2 == 1);
				
    			//std::cout << i << " " << j << " " << k << " : " << g_voxelGrid->isInside(i, j, k) << " and " << (numSurfaceIntersections(coor, direction) % 2 == 1) << "\n";
    		}
    	}
    }
    std::cout << "done trim voxels." << std::endl;
     
    
    //Write out voxel data as obj
    saveVoxelsToObj(argv[2]);
    
    std::cout << "done save voxels obj." << std::endl;
    
    delete g_voxelGrid;
}
