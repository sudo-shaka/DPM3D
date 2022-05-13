#include <cmath>
#include <Cell3D.hpp>
#include <vector>
#include <iostream>
#include <bits/stdc++.h>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <glm/mat3x3.hpp>
#include <glm/gtx/norm.hpp>
#include <thread>
#include <future>

namespace DPM3D{
    //constructors
    Cell::Cell(double _x, double _y, double _z, double _calA0, int f, double _r0 ,double _Kv, double _Ka, double _Kb){
        calA0 = _calA0;
        r0 = _r0;
        Kv = _Kv;
        Ka = _Ka;
        Kb = _Kb;
        int i,j,steps,a,b,c;
        double t = (1+sqrt(5)) / 2;
        Positions.resize(12);
        Positions[0] = glm::dvec3(-1, t, 0);
        Positions[1] = glm::dvec3( 1,  t,  0);
        Positions[2] = glm::dvec3(-1, -t,  0);
        Positions[3] = glm::dvec3( 1, -t,  0);

        Positions[4] = glm::dvec3( 0, -1,  t);
        Positions[5] = glm::dvec3( 0,  1,  t);
        Positions[6] = glm::dvec3( 0, -1, -t);
        Positions[7] = glm::dvec3( 0,  1, -t);

        Positions[8] = glm::dvec3( t,  0, -1);
        Positions[9] = glm::dvec3( t,  0,  1);
        Positions[10] = glm::dvec3(-t,  0, -1);
        Positions[11] = glm::dvec3(-t,  0,  1);
        for(i=0;i<12;i++){
            Positions[i] /= glm::l2Norm(Positions[i]);
        }

        // 5 faces around point 0
        AddFaceIndex(0, 11, 5);
        AddFaceIndex(0, 5, 1);
        AddFaceIndex(0, 1, 7);
        AddFaceIndex(0, 7, 10);
        AddFaceIndex(0, 10, 11);

        // 5 adjacent faces
        AddFaceIndex(1, 5, 9);
        AddFaceIndex(5, 11, 4);
        AddFaceIndex(11, 10, 2);
        AddFaceIndex(10, 7, 6);
        AddFaceIndex(7, 1, 8);

        // 5 faces around point 3
        AddFaceIndex(3, 9, 4);
        AddFaceIndex(3, 4, 2);
        AddFaceIndex(3, 2, 6);
        AddFaceIndex(3, 6, 8);
        AddFaceIndex(3, 8, 9);

        // 5 adjacent faces
        AddFaceIndex(4, 9, 5);
        AddFaceIndex(2, 4, 11);
        AddFaceIndex(6, 2, 10);
        AddFaceIndex(8, 6, 7);
        AddFaceIndex(9, 8, 1);

        std::vector<int> newF;
        newF.resize(3);
        std::vector<std::vector<int>> newFaces;
        for(i=0;i<f;i++){
            steps = FaceIndices.size();
            for(j=0;j<steps;j++){
                a = Cell::AddMiddlePoint(FaceIndices[j][0],FaceIndices[j][1]);
                b = Cell::AddMiddlePoint(FaceIndices[j][1],FaceIndices[j][2]);
                c = Cell::AddMiddlePoint(FaceIndices[j][2],FaceIndices[j][0]);
                newF = newFace(FaceIndices[j][0],a,c);
                newFaces.push_back(newF);
                newF = newFace(FaceIndices[j][1],b,a);
                newFaces.push_back(newF);
                newF = newFace(FaceIndices[j][2],c,b);
                newFaces.push_back(newF);
                newF = newFace(a,b,c);
                newFaces.push_back(newF);
            }
            FaceIndices = newFaces;
            newFaces.clear();
            newFaces.shrink_to_fit();
        }

        NV = (int)Positions.size();
        Velocities.resize(NV);
        Forces.resize(NV);
        ntriangles = FaceIndices.size();

        v0 = (4.0/3.0)*M_PI*pow(r0,3);
        s0 = pow((6*sqrt(M_PI)*v0*calA0),(2.0/3.0));
        a0 = (s0/(double)ntriangles);
        for(i=0;i<NV;i++){
            Positions[i] *= r0;
            Positions[i].x += _x;
            Positions[i].y += _y;
            Positions[i].z += _z;
        }
    }

    void Cell::ResetForces(){
        for(int vi=0;vi<NV;vi++){
            Forces[vi].x = 0.0;
            Forces[vi].y = 0.0;
            Forces[vi].z = 0.0;
            Velocities[vi].x = 0.0;
            Velocities[vi].y = 0.0;
            Velocities[vi].z = 0.0;
        }
    }

    void Cell::ShapeForceUpdate(){
        VolumeForceUpdate();
        AreaForceUpdate();
        BendingForceUpdate();
    }

    void Cell::EulerUpdate(int steps, double dt){
        for(int i=0; i<steps;i++){
            ShapeForceUpdate();
            EulerUpdate(dt);
        }
    }

    void Cell::EulerUpdate(double dt){
        for(int i=0;i<NV;i++){
            Velocities[i] = 0.5*Forces[i];
            Positions[i] += Forces[i]*dt;
        }
        ResetForces();
    }

    void Cell::FIREMinimization(double alpha0, double dt, int itmax, double Ftol){
        int i, ndelay = 20;
        double P, fnorm, fcheck, vnorm, alpha, dtmax,dtmin;
        int npPos, npNeg, fireit;
        P=0;
        fnorm = 0.0;
        vnorm = 0.0;
        alpha = alpha0;
        dtmax = 10.0* dt;
        dtmin = 1e-2*dt;
        double finc = 1.10;
        double fdec = 0.99;
        double falpha = 0.99;
        int nnegmax = 1000;

        npPos = 0.0;
        npNeg = 0.0;

        fireit = 0;
        fcheck = 10*Ftol;

        ResetForces();
        std::cout << "Starting FIRE minimization" << std::endl;
        while((fcheck > Ftol || fireit < ndelay) && fireit < itmax){
            P=0.0;
            for(i=0;i<NV;i++){
                P += glm::dot(Velocities[i],Forces[i]);
            }
            if(P>0){
                npPos++;
                npNeg = 0;
                if (npPos > ndelay && dt*finc < dtmax){
                    dt *= finc;
                }
                alpha *= falpha;
            }
            else{
                npPos = 0;
                npNeg ++;
                if(npNeg > nnegmax){
                    std::cerr << "(!) Error: P > 0 for too long\n";
                    return;
                }

                for(i=0;i<NV;i++){
                    Positions[i] -= 0.5 * dt * Velocities[i];
                    Velocities[i] *= 0.0;
                }

                if(fireit > ndelay && dt * fdec > dtmin){
                    alpha = alpha0;
                    dt *= fdec;
                }
            }

            for(i=0;i<NV;i++){
                Velocities[i] += 0.5 * dt * Forces[i];
            }

            fnorm = 0.0;
            vnorm = 0.0;

            for(i=0;i<NV;i++){
                fnorm += glm::dot(Forces[i],Forces[i]);
                vnorm += glm::dot(Velocities[i],Velocities[i]);
            }
            fnorm = sqrt(fnorm);
            vnorm = sqrt(vnorm);
            if(fnorm > 0){
                for(i=0;i<NV;i++){
                    Velocities[i] = (1-alpha) * Velocities[i] + alpha * (Forces[i]/fnorm)*vnorm;
                }
            }

            for(i=0;i<NV;i++){
                Positions[i] += dt*Velocities[i];
            }

            ShapeForceUpdate();

            for(i=0;i<NV;i++){
                Velocities[i] += 0.5 * dt * Forces[i];
            }
            fcheck = 0.0;
            for(i=0;i<NV;i++){
                fcheck += glm::dot(Forces[i],Forces[i]);
            }
            fcheck = sqrt(fcheck)/NV;
            fireit++;
            if(fireit % 1000 == 0){
                std::cout << "FIRE progress update:" << std::endl;
                std::cout << "	 iterations = " << fireit << std::endl;
                std::cout << "	 fcheck = " << fcheck << std::endl;
                std::cout << "	 dt = " << dt << std::endl;
                std::cout << "	 P = " << P << std::endl;
                std::cout << "	 alpha = " << alpha << std::endl;
            }

        }

        if(fireit == itmax){
            std::cerr << "(!) Fire Minimization did not converge" << std::endl;
            std::cout << "	 iterations = " << fireit << std::endl;
            std::cout << "	 fcheck = " << fcheck << std::endl;
            std::cout << "	 dt = " << dt << std::endl;
            std::cout << "	 P = " << P << std::endl;
            std::cout << "	 alpha = " << alpha << std::endl;
        }
        else{
            std::cout << "FIRE Minimization Converged in "<< fireit << " timesteps" << std::endl;
            std::cout << "	 iterations = " << fireit << std::endl;
            std::cout << "	 fcheck = " << fcheck << std::endl;
            std::cout << "	 dt = " << dt << std::endl;
            std::cout << "	 P = " << P << std::endl;
            std::cout << "	 alpha = " << alpha << std::endl;
        }
    }

    void Cell::VolumeForceUpdate(){
        double volume = GetVolume(),dist, nucdist = (0.5*pow((3*v0)/(4*M_PI),(1/3)));
        double volumeStrain = (volume/v0) - 1.0;
        glm::dvec3 center = GetCOM(),tmp;
        std::vector<int> tri{0,0,0};
        int i,j;
        for(i=0;i<ntriangles;i++){
            tri = FaceIndices[i];
            tmp = glm::cross((Positions[tri[1]]-center)- (Positions[tri[0]]-center),
                (Positions[tri[2]]-center) - (Positions[tri[0]]-center));
            tmp = glm::normalize(tmp);
            for(j=0;j<3;j++){
                dist = distance(center,Positions[tri[j]]);
                if(nucdist < dist)
                    Forces[tri[j]] -= Kv* 0.5 * volumeStrain * (tmp);
                else
                   Forces[tri[j]] -= (1.0-dist/(nucdist))/nucdist * (center-Positions[tri[j]]);
            }
        }
    }
    void Cell::AreaForceUpdate(){
        //double area,areaStrain;
        double length[3],dli[3],dlim1[3],l0 = sqrt((4*a0)/sqrt(3));
        int i,j, im1[] = {2,0,1}, ip1[] = {1,2,0};
        std::vector<int> tri{0,0,0};
        glm::dvec3 center(0.0,0.0,0.0);
        std::vector<glm::dvec3> positions(3,center), lv(3,center), ulv(3,center);
        for(i=0;i<ntriangles;i++){
            tri = FaceIndices[i];
            for(j=0;j<3;j++){
                positions[j] = Positions[tri[j]];
            }
            //center = (positions[0]+positions[1]+positions[2])/3.0;
            for(j=0;j<3;j++){
                lv[j] = positions[ip1[j]] - positions[j];
                length[j] = distance(positions[ip1[j]],positions[j]);
            }
            for(j=0;j<3;j++){
                ulv[j] = lv[j]/length[j];
                dli[j] = length[j]/l0 - 1.0;
                dlim1[j] = (length[im1[j]]/l0) - 1.0;
            }
            for(j=0;j<3;j++){
                Forces[tri[j]] += (Ka*(sqrt(a0)/l0)) * (dli[j]*ulv[j]-dlim1[j]*ulv[im1[j]]);
            }
        }
    }

    void Cell::BendingForceUpdate(){
        int i, j,k,t,c;
        std::vector<std::vector<int>> corner;
        std::vector<int> tri{0,0,0},UsedIndexes,usedPositions;
        glm::dvec3 center;
        double surfaceArea, flatArea,l1,l2,l3,s,bendingStrain;
        for(i=0;i<NV;i++){
            c = 0;
            corner.clear();
            UsedIndexes.clear();
            for(j=0;j<ntriangles;j++){
                tri = FaceIndices[j];
                for(t=0;t<3;t++){
                    if(tri[t] == i){
                        corner.push_back(tri);
                        UsedIndexes.push_back(j);
                        c++;
                        break;
                    }
                }
            }
            surfaceArea = 0.0;
            center *= 0.0;
            usedPositions.clear();

            for(j=0;j<c;j++){
                surfaceArea += GetArea(UsedIndexes[j]);
                for(t=0;t<3;t++){
                    usedPositions.push_back(corner[j][t]);
                }
            }
            sort( usedPositions.begin(), usedPositions.end() );
            usedPositions.erase( unique( usedPositions.begin(), usedPositions.end() ), usedPositions.end() );
            usedPositions.erase(std::remove(usedPositions.begin(), usedPositions.end(), i), usedPositions.end());
            usedPositions.shrink_to_fit();
            k = (int)usedPositions.size();
            for(j=0;j<k;j++){
                center += Positions[usedPositions[j]];
            }
            center /= k;
            flatArea = 0.0;

            for(j=0;j<k;j++){
                if(j != k-1){
                    l1 = distance(center,Positions[usedPositions[j]]);
                    l2 = distance(center,Positions[usedPositions[j+1]]);
                    l3 = distance(Positions[usedPositions[j]],Positions[usedPositions[j+1]]);
                }
                else{
                    l1 = distance(center,Positions[usedPositions[j]]);
                    l2 = distance(center,Positions[usedPositions[0]]);
                    l3 = distance(Positions[usedPositions[j]],Positions[usedPositions[0]]);
                }
                s = (l1+l2+l3)/2;
                flatArea += sqrt(s*(s-l1)*(s-l2)*(s-l3));
            }
            bendingStrain = (surfaceArea/flatArea) -1;
            Forces[i] += Kb*bendingStrain*(glm::normalize(center-Positions[i]));
        }
    }


    std::vector<glm::dvec3> Cell::GetShapeForces(DPM3D::Cell Cell){
        std::vector<glm::dvec3> Forces,Fv,Fa,Fb;
        Forces.resize(Cell.NV);
        Fv.resize(Cell.NV);
        Fa.resize(Cell.NV);
        Fb.resize(Cell.NV);
        auto v = std::async(DPM3D::Cell::GetVolumeForces,Cell); 
        auto a = std::async(DPM3D::Cell::GetAreaForces,Cell);
        auto b = std::async(DPM3D::Cell::GetBendingForces,Cell);
        Fv = v.get();
        Fa = a.get();
        Fb = b.get();
        for(int i=0;i<Cell.NV;i++){
            Forces[i] = Fa[i]+Fb[i]+Fv[i];
        }

        return Forces;
    }
    std::vector<glm::dvec3> Cell::GetVolumeForces(DPM3D::Cell Cell){
        double volume = Cell.GetVolume(),dist, nucdist = (0.5*pow((3*Cell.v0)/(4*M_PI),(1/3)));
        double volumeStrain = (volume/Cell.v0) - 1.0;
        glm::dvec3 center = Cell.GetCOM(),tmp;
        std::vector<glm::dvec3> Forces;
        Forces.resize(Cell.NV);
        std::vector<int> tri{0,0,0};
        int i,j;
        for(i=0;i<Cell.ntriangles;i++){
            tri = Cell.FaceIndices[i];
            tmp = glm::cross((Cell.Positions[tri[1]]-center)- (Cell.Positions[tri[0]]-center),
                (Cell.Positions[tri[2]]-center) - (Cell.Positions[tri[0]]-center));
            tmp = glm::normalize(tmp);
            for(j=0;j<3;j++){
                dist = distance(center,Cell.Positions[tri[j]]);
                if(nucdist < dist)
                    Forces[tri[j]] -= Cell.Kv* 0.5 * volumeStrain * (tmp);
                else
                   Forces[tri[j]] -= (1.0-dist/(nucdist))/nucdist * (center-Cell.Positions[tri[j]]);
            }
        }

        return Forces;
    }
    std::vector<glm::dvec3> Cell::GetAreaForces(DPM3D::Cell Cell){
        double length[3],dli[3],dlim1[3],l0 = sqrt((4*Cell.a0)/sqrt(3));
        int i,j, im1[] = {2,0,1}, ip1[] = {1,2,0};
        std::vector<int> tri{0,0,0};
        glm::dvec3 center(0.0,0.0,0.0);
        std::vector<glm::dvec3> positions(3,center), lv(3,center), ulv(3,center),Forces(Cell.NV,center);
        for(i=0;i<Cell.ntriangles;i++){
            tri = Cell.FaceIndices[i];
            for(j=0;j<3;j++){
                positions[j] = Cell.Positions[tri[j]];
            }
            //center = (positions[0]+positions[1]+positions[2])/3.0;
            for(j=0;j<3;j++){
                lv[j] = positions[ip1[j]] - positions[j];
                length[j] = distance(positions[ip1[j]],positions[j]);
            }
            for(j=0;j<3;j++){
                ulv[j] = lv[j]/length[j];
                dli[j] = length[j]/l0 - 1.0;
                dlim1[j] = (length[im1[j]]/l0) - 1.0;
            }
            for(j=0;j<3;j++){
                Forces[tri[j]] += (Cell.Ka*(sqrt(Cell.a0)/l0)) * (dli[j]*ulv[j]-dlim1[j]*ulv[im1[j]]);
            }
        }
        return Forces;
    }
    std::vector<glm::dvec3> Cell::GetBendingForces(DPM3D::Cell Cell){
        int i, j,k,t,c;
        std::vector<std::vector<int>> corner;
        std::vector<int> tri{0,0,0},UsedIndexes,usedPositions;
        std::vector<glm::dvec3> Forces; Forces.resize(Cell.NV);
        glm::dvec3 center;
        double surfaceArea, flatArea,l1,l2,l3,s,bendingStrain;
        for(i=0;i<Cell.NV;i++){
            c = 0;
            corner.clear();
            UsedIndexes.clear();
            for(j=0;j<Cell.ntriangles;j++){
                tri = Cell.FaceIndices[j];
                for(t=0;t<3;t++){
                    if(tri[t] == i){
                        corner.push_back(tri);
                        UsedIndexes.push_back(j);
                        c++;
                        break;
                    }
                }
            }
            surfaceArea = 0.0;
            center *= 0.0;
            usedPositions.clear();

            for(j=0;j<c;j++){
                surfaceArea += Cell.GetArea(UsedIndexes[j]);
                for(t=0;t<3;t++){
                    usedPositions.push_back(corner[j][t]);
                }
            }
            sort( usedPositions.begin(), usedPositions.end() );
            usedPositions.erase( unique( usedPositions.begin(), usedPositions.end() ), usedPositions.end() );
            usedPositions.erase(std::remove(usedPositions.begin(), usedPositions.end(), i), usedPositions.end());
            usedPositions.shrink_to_fit();
            k = (int)usedPositions.size();
            for(j=0;j<k;j++){
                center += Cell.Positions[usedPositions[j]];
            }
            center /= k;
            flatArea = 0.0;

            for(j=0;j<k;j++){
                if(j != k-1){
                    l1 = distance(center,Cell.Positions[usedPositions[j]]);
                    l2 = distance(center,Cell.Positions[usedPositions[j+1]]);
                    l3 = distance(Cell.Positions[usedPositions[j]],Cell.Positions[usedPositions[j+1]]);
                }
                else{
                    l1 = distance(center,Cell.Positions[usedPositions[j]]);
                    l2 = distance(center,Cell.Positions[usedPositions[0]]);
                    l3 = distance(Cell.Positions[usedPositions[j]],Cell.Positions[usedPositions[0]]);
                }
                s = (l1+l2+l3)/2;
                flatArea += sqrt(s*(s-l1)*(s-l2)*(s-l3));
            }
            bendingStrain = (surfaceArea/flatArea) -1;
            Forces[i] += Cell.Kb*bendingStrain*(glm::normalize(center-Cell.Positions[i]));
        }
        return Forces;
    }

    void Cell::StickToSurface(double z, double mindist){
        glm::dvec3 surfacepos,A,B,com = GetCOM();
        double dist,ftmp;
        std::vector<int> tri{0,0,0};

        for(int i=0;i<ntriangles;i++){
            tri = FaceIndices[i];
            A = Positions[tri[1]] - Positions[tri[0]];
            B = Positions[tri[2]] - Positions[tri[0]];
            for(int j=0; j<3;j++){
                surfacepos = Positions[tri[j]];
                surfacepos.z = z;
                dist = distance(surfacepos,Positions[tri[j]]);
                if((A.x*B.y - A.y*B.x)< 0.0 && dist < mindist){
                    ftmp = (1.0 - dist/(mindist))/mindist;
                    Forces[tri[j]] += ftmp*(glm::normalize(Positions[tri[j]] - com));
                    ExtendVertex(tri[j],a0);
                }
                if(Positions[tri[j]].z < z){
                    Forces[tri[j]].z += 10*pow((Positions[tri[j]].z - z),2);
                }
            }
        }
    }

    void Cell::StickToSurface(double mindist){
        glm::dvec3 A,B,com = GetCOM();
        double dist,ftmp,disttmp;
        std::vector<int> tri{0,0,0};
        int ti,vi,si,siclosest,angle;

        for(ti=0;ti<ntriangles;ti++){
            tri = FaceIndices[ti];
            A = Positions[tri[1]] - Positions[tri[0]];
            B = Positions[tri[2]] - Positions[tri[0]];
            for(vi=0;vi<3;vi++){
                dist = distance(surfacepositions[0],Positions[tri[vi]]);
                siclosest = 0;
                for(si=0;si<nsurfacep;si++){
                    disttmp = distance(surfacepositions[si],Positions[tri[vi]]);
                    if(disttmp < dist){
                        dist = disttmp;
                        siclosest = si;
                    }
                }
                angle = acos(glm::dot(Positions[tri[vi]],com)/(glm::l2Norm(Positions[tri[vi]])*glm::l2Norm(com)));
                if(Positions[tri[vi]].z < surfacepositions[0].z){
                    Forces[tri[vi]].z += 10*pow((Positions[tri[vi]].z - surfacepositions[0].z),2);
                }
                else if((A.x*B.y - A.y*B.x)< 0.0 && dist < mindist){
                    ftmp = (1.0 - dist/(mindist))/mindist;
                    ftmp *= (M_PI/2)-angle;
                    //Forces[tri[vi]] += 0.5*ftmp*(glm::normalize(surfacepositions[siclosest]- Positions[tri[vi]]));
                    Forces[tri[vi]] += 0.5*ftmp*(glm::normalize(surfacepositions[siclosest] - com));
                    //Forces[tri[vi]] += 0.5*ftmp*(glm::normalize(Positions[tri[vi]] - com));
                }
            }
        }
    }

    void Cell::RepelSurface(){
        for(int i =0;i<NV;i++){
            if(Positions[i].z < surfacepositions[0].z){
                Forces[i].z += 10*pow((Positions[i].z - surfacepositions[0].z),2);
            }
        }
    }


    void Cell::ExtendVertex(int vi, double k){
        glm::dvec3 center = GetCOM();
        Forces[vi] += k*glm::normalize(Positions[vi]-center);
    }

    void Cell::SetupSurface(int npoints){
        nsurfacep = npoints;
        glm::dvec3 surfacepos,com = GetCOM();
        int lenp = (int)sqrt(nsurfacep);
        int i,j;
        double maxx=0.0, maxy=0.0;
        double minz=Positions[0].z, miny=Positions[0].y,minx=Positions[0].x;
        double distx,disty,dincx,dincy;
        for(i=0;i<NV;i++){
            if(Positions[i].x < minx){
                minx = Positions[i].x;
            }
            if(Positions[i].y < miny){
                miny = Positions[i].y;
            }
            if(Positions[i].z < minz){
                minz = Positions[i].z;
            }
            if(Positions[i].x > maxx){
                maxx = Positions[i].x;
            }
            if(Positions[i].y > maxy){
                maxy = Positions[i].y;
            }
        }
        distx = (maxx-minx)*5.0;
        disty = (maxy-miny)*5.0;
        dincx = distx/lenp;
        dincy = disty/lenp;
        for(i=0; i<lenp; i++){
            for(j=0;j<lenp;j++){
                surfacepos = glm::dvec3(((-distx/2)+(dincx*(double)i)),((-disty/2)+(dincy*(double)j)),minz);
                surfacepos.x += com.x;
                surfacepos.y += com.y;
                surfacepositions.push_back(surfacepos);
            }
        }
    }

    void Cell::SurfaceStrech(double scale){
        for (int i=0;i<(int)surfacepositions.size();i++){
            surfacepositions[i].x *= scale;
        }
    }

    double Cell::GetVolume(){
        int i,j;
        double det,volume=0.0;
        std::vector<int> tri{0,0,0};
        std::vector<glm::dvec3> positions;
        glm::dmat3 mat;
        glm::dvec3 center = GetCOM();
        for(i=0;i<ntriangles;i++){
            tri = FaceIndices[i];
            positions.clear();
            for(j=0;j<3;j++){
                positions.push_back(Positions[tri[j]]-center);
            }
            mat = glm::dmat3(positions[0],positions[1],positions[2]);
            det = glm::determinant(mat);
            if(det < 0.0)
                det *= -1;
            volume += det;
        }
        volume /= 6;
        return volume;
    }

    double Cell::GetSurfaceArea(){
        double area = 0.0;
        for(int i=0; i<ntriangles;i++){
            area += GetArea(i);
        }
        return area;
    }

    double Cell::GetArea(int i){
        glm::dvec3 a,b,c;
        double lens1,lens2,lens3,s,area;
        a = Positions[FaceIndices[i][0]];
        b = Positions[FaceIndices[i][1]];
        c = Positions[FaceIndices[i][2]];
        lens1 = distance(a,b);
        lens2 = distance(b,c);
        lens3 = distance(c,a);
        s = (lens1+lens2+lens3)/2;
        area = sqrt(s*(s-lens1)*(s-lens2)*(s-lens3));
        return area;
    }

    double Cell::GetCalA0(){
        double s = GetSurfaceArea();
        double v = GetVolume();
        return (pow(s,3/2))/(6*sqrt(M_PI)*v);
    }

    glm::dvec3 Cell::GetCOM(){
        glm::dvec3 compos = glm::dvec3(0.0,0.0,0.0);
        for(int i=0;i<NV;i++){
            compos += Positions[i];
        }
        compos /= NV;
        return compos;
    }

    bool Cell::pointInside(glm::dvec3 point){
        std::vector<int> tri;
        double det, inv_det,u,v;
        glm::dvec3 dir,edge1,edge2,tvec,pvec,qvec;
        int n_crosses = 0;
        dir = glm::normalize(glm::dvec3(1,1,1));
        for(int i=0;i<ntriangles;i++){
            tri = FaceIndices[i];
            edge1 = Positions[tri[1]] - Positions[tri[0]];
            edge2 = Positions[tri[2]] - Positions[tri[0]];
            pvec = glm::cross(dir,edge2);
            det = glm::dot(edge1,pvec);
            /*if(det<ep && det > -ep){
                continue;
            }*/
            inv_det = 1.0/det;
            tvec = point - Positions[tri[0]];
            u = glm::dot(tvec,pvec);
            if(u<0 || u > 1.0 || u*inv_det <0 || u*inv_det > 1.0){
                continue;
            }
            qvec = glm::cross(tvec,edge1);
            v = glm::dot(dir,qvec);
            if(v < 0.0 || u+v > 1.0 || v*inv_det < 0.0 || (u*inv_det)+(v*inv_det) >1.0)
            {
                continue;
            }
            n_crosses++;
        }
        return (n_crosses % 2 != 0);
    }

    //The following functions are for construction and vector maniplulations

    void Cell::AddFaceIndex(int a, int b, int c){
        std::vector<int> idx{a,b,c};
        FaceIndices.push_back(idx);
    }

    int Cell::AddMiddlePoint(int p1, int p2){
        int key; int i;
        if(p1 < p2){
            key = floor((p1+p2) * (p1+p2+1)/2) + p1;
        }
        else{
            key = floor((p1+p2) * (p1+p2+1)/2) + p2;
        }
        for(i=0;i<(int)midpointCache.size();i++){
            if(key == midpointCache[i][0])
                return midpointCache[i][1];
        }

        Cell::AddVertex(GetMiddlePoint(p1,p2));
        i = Positions.size()-1;
        std::vector<int> cache; cache.resize(2);
        cache.shrink_to_fit();
        cache[0] = key; cache[1] = i;
        midpointCache.push_back(cache);
        return i;

    }

    void Cell::AddVertex(glm::dvec3 vec){
        Positions.push_back(vec);
    }


    std::vector<int> newFace(int a, int b, int c){
        std::vector<int> vec{a,b,c};
        return vec;
    }

    glm::dvec3 GetMiddlePoint(glm::dvec3 a, glm::dvec3 b){
        glm::dvec3 pmid = a+b;
        pmid *= 0.5;
        return pmid;
    }

    glm::dvec3 Cell::GetMiddlePoint(int i ,int j){
        glm::dvec3 p1 = Positions[i];
        glm::dvec3 p2 = Positions[j];
        glm::dvec3 pmid = p1+p2;
        pmid  *= 0.5;
        pmid /= glm::l2Norm(pmid);
        return pmid;
    }

    double distance(glm::dvec3 a,glm::dvec3 b){
        return sqrt(distanceSq(a,b));
    }
    double distanceSq(glm::dvec3 a,glm::dvec3 b){
        glm::dvec3 temp = a-b;
        return glm::dot(temp,temp);
    }

}
