/*! \file main.cxx
 *  \brief Main program

*/

#include "TreeFrog.h"

using namespace std;
using namespace Math;
using namespace NBody;

int main(int argc,char **argv)
{
#ifdef USEMPI
    //start MPI
    MPI_Init(&argc,&argv);
    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&NProcs);
    //and this processes' rank is
    MPI_Comm_rank(MPI_COMM_WORLD,&ThisTask);
#else
    int ThisTask=0,NProcs=1,Nlocal,Ntotal;
    int StartSnap,EndSnap;
#endif
    int nthreads;
//#ifdef USEOPENMP
//#pragma omp parallel
//    {
//    if (omp_get_thread_num()==0) nthreads=omp_get_num_threads();
//    }
//#else
    nthreads=1;
//#endif

#ifdef USEMPI
    if (ThisTask==0) cout<<"VELOCIraptor/Tree running with MPI. Number of mpi threads: "<<NProcs<<endl;
#endif
//#ifdef USEOPENMP
//    if (ThisTask==0) cout<<"VELOCIraptor/Tree running with OpenMP. Number of openmp threads: "<<nthreads<<endl;
//#endif

    //store the configuration settings in the options object
    Options opt;

    //store the halo tree
    HaloTreeData *pht;
    //store the progenitors of a given set of objects
    ProgenitorData **pprogen;
    //if multiple snapshots are used in constructing the progenitor list, must store the possible set of descendents
    //which can then be cleaned to ensure one descendent
    ///\todo eventually must clean up graph construction so that the graph is completely reversible
    DescendantDataProgenBased **pprogendescen;

    //stoe the descendants of a given set of objects
    DescendantData **pdescen;
    //if multiple snapshots are used in constructing the descendant list, must store the possible set of progenitors
    //which can then be cleaned to ensure one progenitor
    ProgenitorDataDescenBased **pdescenprogen;
    ProgenitorDataDescenBased *pdescenprogentemp;
    //these are used to store temporary sets of posssible candidate progenitors/descendants when code links across multiple snaps
    ProgenitorData *pprogentemp;
    DescendantData *pdescentemp;
    //store the halo ids of particles in objects, with mapping of ids, makes it easy to look up the haloid of a given particle
    unsigned int *pfofp,*pfofd;
    unsigned int *prank=NULL;
    //int *ibuildrankflag;
    long long i,j;
    long unsigned nh,nhp,nhd;
    long unsigned newnp;
    //flags indicating whether updates to list are necessary
    int ilistupdated;

    //store the time taken by segments of the code
    double time1,time2;

    GetArgs(argc, argv, opt);
#ifdef USEMPI
    MPILoadBalanceSnapshots(opt);
#else
    StartSnap=0;
    EndSnap=opt.numsnapshots;
#endif
    //read particle information and allocate memory
    if (ThisTask==0) cout<<"Read Data ... "<<endl;
    pht=ReadData(opt);
    if (ThisTask==0) {
        cout<<"Done Loading"<<endl;
        cout<<"Found "<<opt.TotalNumberofHalos<<" halos "<<endl;
    }
    if (opt.iverbose) {
        long unsigned sum=0;
        for (i=StartSnap;i<EndSnap;i++)
            for (j=0;j<pht[i].numhalos;j++) sum+=pht[i].Halo[j].NumberofParticles;
        cout<<ThisTask<<" has allocated at least "<<sum*sizeof(IDTYPE)/1024./1024./1024.<<"GB of memory to store particle ids"<<endl;
    }

    //computationally expensive but memory efficient mapping of data
    if (opt.imapping==DMEMEFFICIENTMAP) {
        MemoryEfficientMap(opt,pht);
    }
    else {
        //adjust ids if particle ids need to be mapped to index
        if (opt.imapping>DNOMAP) MapPIDStoIndex(opt,pht);
        //check that ids are within allowed range for linking
        IDcheck(opt,pht);
    }
    if (ThisTask==0) {
        cout<<"Memory needed to store addressing for ids "<<sizeof(IDTYPE)*opt.MaxIDValue/1024./1024.0/1024.0<<"GB, maximum ID of "<<opt.MaxIDValue<<endl;
        if (opt.imerittype==MERITRankWeightedBoth) cout<<" And also need similar amount of memory to store ranking information needed by merit"<<endl;
    }

// ASSUME A DESCENDANT BASED TREE IS USED
    time1=MyGetTime();

    if (opt.iverbose) cout<<"Starting descendant cross matching "<<endl;

    pdescen=new DescendantData*[opt.numsnapshots];
    pfofd=new unsigned int[opt.MaxIDValue];
    
    for (i=0;i<opt.MaxIDValue;i++) pfofd[i]=0;
    
    if (opt.imerittype==MERITRankWeightedBoth) {
        prank=new unsigned int[opt.MaxIDValue];
        for (i=0;i<opt.MaxIDValue;i++) prank[i]=0;
    }
    
    //Store descendant based progenitor information
    pdescenprogen=new ProgenitorDataDescenBased*[opt.numsnapshots];
    for (i=0;i<opt.numsnapshots;i++) pdescenprogen[i]=NULL;

/* SAVE FIRST STEP */
    if (pht[StartSnap].numhalos>0) pdescenprogen[StartSnap]=new ProgenitorDataDescenBased[pht[StartSnap].numhalos];
    
 /* FIRST DO SINGLE STEP LINKING */
    for (i=0;i<opt.numsnapshots;i++) {
        if (!(i>=StartSnap && i<EndSnap-1)) {pdescen[i]=NULL; continue;}
        time2=MyGetTime();
        if(pht[i].numhalos==0) {pdescen[i]=NULL;continue;}
        cout<<i<<" "<<pht[i].numhalos<<" cross matching objects in descendant direction"<<endl;
        
        for (j=1;j<=opt.numsteps;j++)
            if (i+j<EndSnap)
                if (pdescenprogen[i+j]==NULL && pht[i+j].numhalos>0)
                    pdescenprogen[i+j]=new ProgenitorDataDescenBased[pht[i+j].numhalos];
        
        /* if beyond endsnap allocate mem but do nothing */
        if (i>=EndSnap) {
            pdescen[i]=new DescendantData[pht[i].numhalos];
            continue;
        }
        
        Int_t istep=1;
        if (!(i+istep<=opt.numsnapshots-1)) continue;
        if (!(i+istep<EndSnap)) continue;
        
        /* Set pfof progenitor data structure, used to produce links. Only
        produced IF snapshot not first one */

        for (j=0;j<pht[i+istep].numhalos;j++) {
            for (int k=0;k<pht[i+istep].Halo[j].NumberofParticles;k++) {
                pfofd[pht[i+istep].Halo[j].ParticleID[k]]=j+1;
            }
        }

        if (opt.imerittype==MERITRankWeightedBoth) {
            for (j=0;j<pht[i+istep].numhalos;j++) {
                for (int k=0;k<pht[i+istep].Halo[j].NumberofParticles;k++) {
                    prank[pht[i+istep].Halo[j].ParticleID[k]]=k+1;
                }
            }
        }
            
        /* Now if also doing core weighting then update the halo id associated
         with the particle so that it is its current halo core ID + total number
         of halos */
        if (opt.icorematchtype!=PARTLISTNOCORE && opt.particle_frac<1 && opt.particle_frac>0) {
            for (j=0;j<pht[i+istep].numhalos;j++) {
                newnp=max((Int_t)(pht[i+istep].Halo[j].NumberofParticles*opt.particle_frac),opt.min_numpart);
                newnp=min(pht[i+istep].Halo[j].NumberofParticles,newnp);
                for (Int_t k=0;k<newnp;k++)
                    pfofd[pht[i+istep].Halo[j].ParticleID[k]]=j+1+pht[i+istep].numhalos;
            }
        }

        /* Begin cross matching with snapshot(s) for first linking, cross match and
         allocate memory identify candidate descendants */
        pdescen[i]=CrossMatchDescendant(opt,  pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pfofd, ilistupdated, istep, prank);
        /* Update halo ids */
        UpdateDescendantIndexing(istep, pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pdescen[i]);
        /* Build a temporally local descendant based progenitor data */
        BuildDescendantBasedProgenitorList(i, pht[i].numhalos, pdescen[i], pdescenprogen[i+istep]);
        /* Rank progenitors at time i of descedants found at time i+istep based on
         their merit. Ranking is necessary to determine main/secondary branches */
        UpdateDescendantUsingDescendantBasedProgenitorList(pht[i+istep].numhalos, pdescen[i], pdescenprogen[i+istep], istep, opt.meritlimit);
        /* Clean up the information stored in this list, adjusing rankings as necessary */
        CleanCrossMatchDescendant(opt, i, pht, pdescenprogen, pdescen);

        for (j=0;j<pht[i+istep].numhalos;j++) {
            for (int k=0;k<pht[i+istep].Halo[j].NumberofParticles;k++) {
                pfofd[pht[i+istep].Halo[j].ParticleID[k]]=0;
            }
        }

        if (opt.numsteps==1) {   /* Free memory - not required in this case*/
            for (j=0;j<pht[i].numhalos;j++) {
                delete[] pht[i].Halo[j].ParticleID;pht[i].Halo[j].ParticleID=NULL;
            }
        }

        if (opt.iverbose) cout<<ThisTask<<" Finished first pass for descendant processing for snapshot "<<i<<" in "<<MyGetTime()-time2<<endl;

    } // for (i=0;i<opt.numsnapshots;i++) {
// END OF SINGLE STEP LINKING

        
//        // ... THEN LINK OVER MULTIPLE SNAPSHOTS
//        if (opt.numsteps>1) {
//
//                //followed by multi-snapshot linking

//            for (i=0;i<opt.numsnapshots;i++) {
//                if (!(i>=StartSnap && i<EndSnap-1)) continue;
//                time2=MyGetTime();
//                if(pht[i].numhalos==0) continue;
//                cout<<i<<" "<<pht[i].numhalos<<" cross matching objects in descendant direction second pass for poor/missing matches "<<endl;
//                for (j=2;j<=opt.numsteps;j++) if (i+j<EndSnap)
//                    if (pdescenprogen[i+j]==NULL && pht[i+j].numhalos>0) pdescenprogen[i+j]=new ProgenitorDataDescenBased[pht[i+j].numhalos];
//                if (i>=EndSnap) continue;

//                for (Int_t istep=2;istep<=opt.numsteps;istep++) {
//                    if (opt.imerittype==MERITRankWeightedBoth) {
//                        for (j=0;j<pht[i+istep].numhalos;j++) {
//                            for (int k=0;k<pht[i+istep].Halo[j].NumberofParticles;k++) {
//                                prank[pht[i+istep].Halo[j].ParticleID[k]]=k+1;
//                            }
//                        }
//                    }
//                    //now if also doing core weighting then update the halo id associated with the particle so that
//                    //it is its current halo core ID + total number of halos
//                    if (opt.icorematchtype!=PARTLISTNOCORE && opt.particle_frac<1 && opt.particle_frac>0) {
//                        for (j=0;j<pht[i+istep].numhalos;j++) {
//                            newnp=max((Int_t)(pht[i+istep].Halo[j].NumberofParticles*opt.particle_frac),opt.min_numpart);
//                            newnp=min(pht[i+istep].Halo[j].NumberofParticles,newnp);
//                            for (Int_t k=0;k<newnp;k++)
//                                pfofd[pht[i+istep].Halo[j].ParticleID[k]]=j+1+pht[i+istep].numhalos;
//                        }
//                    }
//
//                    //if more than a single step is used to find descendants then we first search i+istep but only for those haloes that are deemed to have
//                    //less than ideal descendants.
//                    pdescentemp=CrossMatchDescendant(opt, pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pfofd, ilistupdated, istep, prank, pdescen[i]);
//                    //if some new descendants are found then need to clean-up and merge information
//                    if (ilistupdated>0) {
//                        //update the halo ids
//                        UpdateDescendantIndexing(istep, pht[i].numhalos, pht[i+istep].numhalos, pht[i].Halo, pht[i+istep].Halo, pdescentemp);
//                        /*
//                        //to rank progenitors of descendants at this time, need to allocate a ProgenitorDataDescenBased list
//                        pdescenprogentemp=new ProgenitorDataDescenBased[pht[i+istep].numhalos];
//                        BuildDescendantBasedProgenitorList(i, pht[i].numhalos, pdescentemp, pdescenprogentemp, istep);
//                        UpdateDescendantUsingDescendantBasedProgenitorList(pht[i+istep].numhalos, pdescentemp, pdescenprogentemp, istep, opt.meritlimit);
//                        //having ranked the progenitors based on their descendants looking backwards, we can now update the descendant list appropriately
//                        UpdateRefDescendants(opt,pht[i].numhalos, pdescen[i], pdescentemp, pdescenprogen, i);
//                        //clean up the information stored in this list, adjusing rankings as necessary
//                        CleanCrossMatchDescendant(opt, i, pht, pdescenprogen, pdescen);
//                        delete[] pdescenprogentemp;
//                        */
//                        BuildDescendantBasedProgenitorList(i, pht[i].numhalos, pdescentemp, pdescenprogen[i+istep], istep);
//                        UpdateDescendantAcrossTimeUsingDescendantBasedProgenitorList(pht[i+istep].numhalos, pdescen, pdescentemp, pdescenprogen[i+istep], istep, opt.meritlimit);
//                        //having ranked the progenitors based on their descendants looking backwards, we can now update the descendant list appropriately
//                        UpdateRefDescendants(opt,pht[i].numhalos, pdescen[i], pdescentemp, pdescenprogen, i);
//                        //clean up the information stored in this list, adjusing rankings as necessary
//                        CleanCrossMatchDescendant(opt, i, pht, pdescenprogen, pdescen);
//                    }
//                    delete[] pdescentemp;
//
//                    for (j=0;j<pht[i+istep].numhalos;j++) {
//                        for (int k=0;k<pht[i+istep].Halo[j].NumberofParticles;k++) {
//                            pfofd[pht[i+istep].Halo[j].ParticleID[k]]=0;
//                        }
//                    }
//                } // for (Int_t istep=2;istep<=opt.numsteps;istep++) {
//                //to free up some memory, no need to keep particle ids
//                for (j=0;j<pht[i].numhalos;j++) {
//                    delete[] pht[i].Halo[j].ParticleID;pht[i].Halo[j].ParticleID=NULL;
//                }
//                if (opt.iverbose) cout<<ThisTask<<" finished descendant processing for snapshot "<<i<<" in "<<MyGetTime()-time2<<endl;
//            } // for (i=0;i<opt.numsnapshots;i++) {
            /*
            //and final reranking where descendant looks at all progenitors
            for (i=0;i<opt.numsnapshots;i++) {
                if (!(i>=StartSnap+1 && i<EndSnap))  continue;
                UpdateDescendantAcrossTimeUsingDescendantBasedProgenitorList(pht[i].numhalos, pdescen, pdescenprogen[i], opt.meritlimit);
            }
            */
//        } //if (opt.numsteps>1) {

//        delete[] pfofd;
//        if (opt.imerittype==MERITRankWeightedBoth) delete[] prank;
//
//        //after first pass (and localized to mpi domains) start final cleaning of descendant tree.
//#ifdef USEMPI
//        //note that for full cleanning the simpliest mpi overlap ???
//        if (NProcs>1) {
//            MPIUpdateDescendantUsingProgenitors(opt, pht, pdescenprogen);
//            if (opt.numsteps>1) MPIUpdateDescendants(opt, pht, pdescen);
//
//        }
//#endif
//        //first we rank the progenitors in the descendant tree using the descedants
//        for (i=0;i<opt.numsnapshots;i++) {
//            //check if data is load making sure i is in appropriate range (note that only look above StartSnap (as first snap can't have progenitors)
//            //not that last snapshot locally has halos with no descendants so no need to clean this list
//            if (i>=StartSnap && i<EndSnap) {
//                if (opt.iverbose) cout<<ThisTask<<" Rank progenitors in descendant tree using descendants "<<i<<endl;
//                if (pdescenprogen[i]!=NULL) {
//                    RankDescendantProgenitors(i, pht, pdescenprogen, pdescen, opt.iopttemporalmerittype);
//                }
//            }
//        }
//        //clean the descendant tree of progenitors with multiple primary ranked descendants
//        for (i=0;i<opt.numsnapshots-1;i++) {
//            if (i>=StartSnap && i<EndSnap-1) {
//                if (opt.iverbose) cout<<ThisTask<<" Cleaning descendant tree of progenitors that are considered primary progenitors of multiple descendants "<<i<<endl;
//                CleanCrossMatchDescendant(opt, i, pht, pdescenprogen, pdescen);
//            }
//        }
//        //then clean descendant tree for any objects with no primary ranked progenitors
//        for (i=0;i<opt.numsnapshots;i++) {
//            if (i>=StartSnap+1 && i<EndSnap) {
//                if (opt.iverbose) cout<<ThisTask<<" Cleaning descendant tree for missing progenitors "<<i<<endl;
//                if (pdescenprogen[i]!=NULL) {
//                    CleanDescendantsForMissingProgenitors(opt, i, pht, pdescenprogen, pdescen);
//                }
//            }
//        }
//        //free some memory
//        for (i=0;i<opt.numsnapshots;i++) {
//            if (i>=StartSnap && i<EndSnap)
//            {
//                delete[] pdescenprogen[i];
//                pdescenprogen[i]=NULL;
//            }
//        }
//        if (opt.iverbose) cout<<"Finished Descendant cross matching "<<MyGetTime()-time1<<endl;
//        //finally rerank the descendant list based on the ranking stored in the desecndant to progenitor value and then merit.
//        time1=MyGetTime();
//        RerankDescendants(opt,pht,pdescen);
//        if (opt.iverbose) cout<<"Finished reranking "<<MyGetTime()-time1<<endl;

//now adjust the halo ids as prior, halo ids simple represent read order, allowing for addressing of memory by halo id value
//    UpdateHaloIDs(opt,pht);
#ifdef USEMPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

/* WRITE RESULTS TO FILE */
    if (opt.icatalog!=DCROSSCAT) {
        if (opt.icatalog!=DGRAPH) {
            WriteHaloMergerTree(opt,pdescen,pht);
        } else {
            char fname[1000];
            sprintf(fname,"%s.graph",opt.outname);
            sprintf(opt.outname,"%s",fname);
            WriteHaloGraph(opt,pprogen,pdescen,pht);
        }
    }
    else WriteCrossComp(opt,pprogen,pht);


        
// FREE MEMORY
    if (opt.iverbose) cout<<"Freeing memory..."<<endl;

    for (i=0;i<opt.numsnapshots;i++)
    {
        if (pdescen[i]!=NULL)
            delete[] pdescen[i];
    }
    delete[] pdescen;

    if (opt.numsteps>1) {
        for (i=0;i<opt.numsnapshots;i++) {
    #ifdef USEMPI
            //check if data is load making sure i is in appropriate range
                if (i>=StartSnap && i<EndSnap) {
    #endif
            if (pdescenprogen[i]!=NULL)
                delete[] pdescenprogen[i];
    #ifdef USEMPI
            }
   #endif
            }
        }
    delete[] pdescenprogen;
    
// Clean up MPI, if used
    
#ifdef USEMPI
    delete[] mpi_startsnap;
    delete[] mpi_endsnap;
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    delete[] pht;

#ifdef USEMPI
    MPI_Finalize();
#endif
    return 0;
}
