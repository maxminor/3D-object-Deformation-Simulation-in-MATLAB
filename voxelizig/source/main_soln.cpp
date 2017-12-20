//Computational Fabrication Assignment #1
#include <iostream>
#include <vector>
#include "../include/CompFab.h"
#include "../include/Mesh.h"

//Ray-Triangle Intersection
//Returns 1 if triangle and ray intersect, 0 otherwise
int rayTriangleIntersection(CompFab::Ray &ray, CompFab::Triangle &triangle)
{
    CompFab::Vec3 e1 = triangle.m_v2 - triangle.m_v1;
    CompFab::Vec3 e2 = triangle.m_v3 - triangle.m_v1;
    CompFab::Vec3 P, Q, T;
    float det, inv_det, u, v;
    float t;

    //Begin calculating determinant - also used to calculate u parameter
    P =  ray.m_direction%e2;

    //if determinant is near zero, ray lies in plane of triangle
    det = e1*P;

    //NOT CULLING
    if(det > -EPSILON && det < EPSILON)
        return 0;

    inv_det = 1.f / det;

    //calculate distance from V1 to ray origin
    T = ray.m_origin - triangle.m_v1;

    //Calculate u parameter and test bound
    u = T*P*inv_det;

    //The intersection lies outside of the triangle
    if(u < 0.f || u > 1.f) return 0;

    //Prepare to test v parameter
    Q = T%e1;

    //Calculate V parameter and test bound
    v = ray.m_direction*Q*inv_det;

    //The intersection lies outside of the triangle
    if(v < 0.f || u + v  > 1.f)
        return 0;

    t = e2*Q*inv_det;

    if(t > EPSILON) { //ray intersection
        return 1;
    }

    // No hit, no win
    return 0;

}

//Triangle list (global)
typedef std::vector<CompFab::Triangle> TriangleList;

TriangleList g_triangleList;
CompFab::VoxelGrid *g_voxelGrid;

//Number of intersections with surface made by a ray originating at voxel and cast in direction.
int numSurfaceIntersections(CompFab::Vec3 &voxel, CompFab::Vec3 &direction)
{

    unsigned int numHits = 0;

    //Cast ray in x-direction
    CompFab::Ray ray;
    ray.m_origin = voxel;
    ray.m_direction = direction;

    //Intersect ray with all triangles and count (NOTE: This is very silly)
    for(unsigned int tri=0; tri<g_triangleList.size(); ++tri)
    {
        numHits += rayTriangleIntersection(ray, g_triangleList[tri]);
    }

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

    unsigned int dim = 64; //dimension of voxel grid (e.g. 32x32x32)

    //Load OBJ
    if(argc < 3)
    {
        std::cout<<"Usage: Voxelizer InputMeshFilename OutputMeshFilename \n";
        //exit(0);
    }

    std::cout<<"Load Mesh : "<<argv[1]<<"\n";
    loadMesh(argv[1], dim);


    std::cout<<"Voxelizing with dimension "<<dim<<"\n";
    //Cast ray, check if voxel is inside or outside
    //even number of surface intersections = outside (OUT then IN then OUT)
    // odd number = inside (IN then OUT)
    CompFab::Vec3 voxelPos;
    CompFab::Vec3 direction(1.0,0.0,0.0);

    for(unsigned int voxelX = 0; voxelX < g_voxelGrid->m_dimX; ++voxelX)
            for(unsigned int voxelY = 0; voxelY < g_voxelGrid->m_dimY; ++voxelY)
                    for(unsigned int voxelZ = 0; voxelZ < g_voxelGrid->m_dimZ; ++voxelZ)
                    {
                        voxelPos[0] = g_voxelGrid->m_lowerLeft[0] + g_voxelGrid->m_spacing*voxelX;
                        voxelPos[1] = g_voxelGrid->m_lowerLeft[1] + g_voxelGrid->m_spacing*voxelY;
                        voxelPos[2] = g_voxelGrid->m_lowerLeft[2] + g_voxelGrid->m_spacing*voxelZ;

                        if((numSurfaceIntersections(voxelPos, direction) % 2) == 0)
                        {
                            g_voxelGrid->isInside(voxelX, voxelY, voxelZ) = false;
                        } else {
                            g_voxelGrid->isInside(voxelX, voxelY, voxelZ) = true;
                        }

                    }

    //Write out voxel data as obj
    saveVoxelsToObj(argv[2]);

    delete g_voxelGrid;
}
