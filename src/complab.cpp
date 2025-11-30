/* This file is a part of the CompLaB3D program.
 *
 * Original 2D CompLaB developed since 2022 by:
 *   University of Georgia (United States)
 *   Chungnam National University (South Korea)
 *
 * 2D Version Developer:
 *   Heewon Jung
 *   Department of Geological Sciences, Chungnam National University
 *   hjung@cnu.ac.kr
 *
 * 3D Version (CompLaB3D) Developers:
 *   Shahram Asgari (shahram.asgari@uga.edu)
 *   Christof Meile (cmeile@uga.edu)
 *   Department of Marine Sciences, University of Georgia
 *
 * https://github.com/shahram444/CompLaB3D-Code
 *
 * CompLaB3D is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "complab_functions.hh"
#include "complab_processors.hh"

#include <chrono>
#include <string>
#include <iostream>
#include <cstring>
#include <vector>
#include <sys/stat.h>

typedef double T;

int main(int argc, char** argv) {
    plbInit(&argc, &argv);
    global::timer("total").start();
    ImageWriter<T> image("leeloo");

    plint kns_count = 0, lb_count = 0;
    char* main_path = (char*)malloc(100 * sizeof(char));
    getcwd(main_path, 100 * sizeof(char));
    char* src_path = (char*)malloc(100 * sizeof(char));
    char* input_path = (char*)malloc(100 * sizeof(char));
    char* output_path = (char*)malloc(100 * sizeof(char));
    char* ns_filename = (char*)malloc(100 * sizeof(char));
    plint nx, ny, nz, num_of_substrates;
    T dx, dy, dz, deltaP, Pe, charcs_length;
    std::string geom_filename, mask_filename;
    std::vector<bool> vec_left_btype, vec_right_btype;
    std::vector<T> vec_c0, vec_left_bcondition, vec_right_bcondition;
    std::vector< std::vector<int>> vec_const_loc;
    std::vector< std::vector<T>>  vec_Kc, vec_Kc_kns, vec_const_lb, vec_const_ub;

    std::string ade_filename;

    bool read_NS_file = 0, read_ADE_file = 0, soluteDindex = 0, track_performance = 0., halfflag = 0, eqflag = 0;

    plint no_dynamics = 0, bounce_back = 1, ns_rerun_iT0 = 0, ns_update_interval = 1, ade_update_interval = 1,
    ns_maxiTer_1, ns_maxiTer_2, ade_rerun_iT0 = 0, ade_maxiTer = 10000000, ade_VTI_iTer = 1000, ade_CHK_iTer = 1000000;

    T tau = 0.6, ns_converge_iT1 = 1e-8, ns_converge_iT2 = 1e-4, ade_converge_iT = 1e-8;

    std::vector<bool> vec_fixC, vec_fixLB;
    std::vector<plint> pore_dynamics, solver_type, reaction_type, vec_sense;

    std::vector<T> vec_solute_poreD;
    std::vector<std::string> vec_subs_names;

    std::string str_mainDir = main_path;

    if (std::to_string(str_mainDir.back()).compare("/") != 0) { str_mainDir += "/"; }

    std::srand(std::time(nullptr));

//-------------------- load the input file and initialize all input parameters --------------------//

    int erck = 0;
    try {
        erck = initialize_complab(
        main_path, src_path, input_path, output_path, ns_filename, ade_filename, geom_filename, mask_filename,
        read_NS_file, ns_rerun_iT0, ns_converge_iT1, ns_converge_iT2, ns_maxiTer_1, ns_maxiTer_2, ns_update_interval, ade_update_interval,
        read_ADE_file, ade_rerun_iT0, ade_VTI_iTer, ade_CHK_iTer, ade_converge_iT, ade_maxiTer, nx, ny, nz, dx, dy, dz, deltaP, tau,
        Pe, charcs_length, vec_solute_poreD, soluteDindex, pore_dynamics, bounce_back, no_dynamics, num_of_substrates, vec_subs_names,
        solver_type, lb_count, kns_count, reaction_type, vec_c0, vec_left_btype, vec_right_btype, vec_left_bcondition, vec_right_bcondition,
        vec_Kc, vec_Kc_kns, vec_fixLB, vec_fixC, vec_sense, vec_const_loc, vec_const_lb, vec_const_ub, track_performance, halfflag, eqflag);

    }

    catch (PlbIOException& exception) {
        pcout << exception.what() << " Terminating the simulation.\n" << std::endl;
        return -1;
    }

    if (erck != 0) { 
        return -1; 
    }

    plint rxn_count = kns_count;

    std::string  str_inputDir = input_path, str_outputDir = output_path;

    if (std::to_string(str_inputDir.back()).compare("/") != 0) { str_inputDir += "/"; }

    if (std::to_string(str_outputDir.back()).compare("/") != 0) { str_outputDir += "/"; }


    /*  =================================== NS Lattice setup ===================================  */

    struct stat statStruct;

    stat(output_path, &statStruct);

    pcout << "CompLaB3D main directory = " << str_mainDir << std::endl;
    pcout << "CompLaB3D input directory = " << main_path << "/" << input_path << std::endl;
    pcout << "CompLaB3D output directory = " << main_path << "/" << output_path << std::endl << std::endl;

    if (S_ISDIR(statStruct.st_mode)) {}
    else { mkdir(output_path, 0777); }

    global::directories().setOutputDir(str_outputDir);

    T PoreMeanU = 0, PoreMaxUx = 0, DarcyOutletUx = 0, DarcyMiddleUx =0, DarcyInletUx =0;

    plint iT = 0;

    T nsLatticeTau = tau;

    T nsLatticeOmega = 1 / nsLatticeTau;

    T nsLatticeNu = DESCRIPTOR<T>::cs2 * (nsLatticeTau - 0.5);
    T nsNu = DESCRIPTOR<T>::cs2 * (nsLatticeTau - 0.5);
    pcout << "Kinematic Viscosity (nsNu) = " << nsNu << ".\n";

    char* ns_read_filename = strcat(strdup(str_inputDir.c_str()), ns_filename);

    MultiScalarField3D<int> geometry(nx, ny, nz);

    readGeometry(str_inputDir + geom_filename, geometry);

    saveGeometry("inputGeom", geometry);

    MultiScalarField3D<int> distanceDomain(nx, ny, nz);

    distanceDomain = geometry;

    std::vector< std::vector<std::vector<plint> > > distVec(nx, std::vector<std::vector<plint>>(ny, std::vector<plint>(nz)));

    calculateDistanceFromSolid3D(distanceDomain, no_dynamics, bounce_back, distVec);

    applyProcessingFunctional(new createDistanceDomain3D<int>(distVec), distanceDomain.getBoundingBox(), distanceDomain);

    MultiScalarField3D<int> ageDomain(nx, ny, nz);

    ageDomain = geometry;

    applyProcessingFunctional(new createAgeDomain3D<int>(pore_dynamics, bounce_back, no_dynamics), ageDomain.getBoundingBox(), ageDomain);

    if (track_performance == 1) {
        pcout << "Performance tracker has been activated. Skipping all the non-essential IO.\n";
    }

    pcout << "\nInitializing the fluid lattice (deltaP = " << deltaP << ").\n";

    MultiBlockLattice3D<T, DESCRIPTOR> nsLattice(nx, ny, nz, new IncBGKdynamics<T, DESCRIPTOR>(nsLatticeOmega));

    util::ValueTracer<T> ns_convg1(1.0, 1000.0, ns_converge_iT1);

    NSdomainSetup(nsLattice, createLocalBoundaryCondition3D<T, DESCRIPTOR>(), geometry, deltaP, nsLatticeOmega, pore_dynamics, bounce_back, no_dynamics);

    //-------------------- NS Lattice main loop --------------------//

    global::timer("NS").start();

        if (Pe == 0) {
            pcout << "Peclet number is set to 0. Skipping a lattice Boltzmann flow solver.\n";
        }

        else {
            pcout << "nsLatticeTau = " << nsLatticeTau << ", nsLatticeOmega = " << nsLatticeOmega << ", nsLatticeNu = " << nsLatticeNu << std::endl;
            pcout << "\n========== LBM NS simulation begins ==========\n";

            if (read_NS_file == 1 && track_performance == 0) {
                pcout << "run continuous simulation for nsLattice." << std::endl;
                pcout << "Load binary block for nsLattice." << std::endl;

                try { 
                    loadBinaryBlock(nsLattice, strcat(ns_read_filename, ".chk")); 
                }
                catch (PlbIOException& exception) { 
                    pcout << exception.what() << ". Terminating the simulation.\n" << std::endl; 
                    return -1; 
                }

                if (ns_rerun_iT0 == 0) { 
                    pcout << "Use the existing checkpoint file for nsLattice." << std::endl; 
                }
                else {
                    pcout << "run the main nsLattice loop to a new steady state" << std::endl;
                    
                    iT = ns_rerun_iT0;
                    for (; iT < ns_maxiTer_1; ++iT) {
                        nsLattice.collideAndStream();
                        ns_convg1.takeValue(getStoredAverageEnergy(nsLattice), true);
                        if (ns_convg1.hasConverged()) { break; }
                    }
                }

            }

            else {
                pcout << "Run a new simulation for nsLattice" << std::endl;

                for (; iT < ns_maxiTer_1; ++iT) {
                    nsLattice.collideAndStream();
                    ns_convg1.takeValue(getStoredAverageEnergy(nsLattice), true);
                    if (ns_convg1.hasConverged()) { break; }
                }
            }
        }
      
        PoreMeanU = computeAverage( *computeVelocityNorm(nsLattice, Box3D (1,nx-2,0,ny-1,0,nz-1)) );
        pcout << "PoreMeanU = " << PoreMeanU << ".\n";

        PoreMaxUx = computeMax(*computeVelocityComponent(nsLattice, Box3D(1, nx - 2, 0, ny - 1, 0, nz - 1), 0));
        pcout << "PoreMaxUx = " << PoreMaxUx << ".\n";

        DarcyOutletUx = computeAverage(*computeVelocityComponent(nsLattice, Box3D(nx - 2, nx - 2, 0, ny - 1, 0, nz - 1), 0));
        pcout << "DarcyOutletUx = " << DarcyOutletUx << ".\n";

        DarcyMiddleUx = computeAverage(*computeVelocityComponent(nsLattice, Box3D((nx - 1) / 2, (nx - 1) / 2, 0, ny - 1, 0, nz - 1), 0));
        pcout << "DarcyMiddleUx = " << DarcyMiddleUx << ".\n";

        DarcyInletUx = computeAverage(*computeVelocityComponent(nsLattice, Box3D(1, 1, 0, ny - 1, 0, nz - 1), 0));
        pcout << "DarcyInletUx = " << DarcyInletUx << ".\n";

        T Ma = PoreMaxUx / sqrt(BGK<T>::cs2);

        pcout << "CFL number (= maximum local lattice velocity)= " << PoreMaxUx << ".\n";

        pcout << "Mach number = " << Ma << ".\n";

        if (Ma > 1) { pcout << "Ma must be << 1. Terminating the simulation." << std::endl; return -1; }        

        // Calculate initial permeability
        T permeability = computePermeability(nsLattice, nsNu, deltaP, nsLattice.getBoundingBox());
        pcout << "Initial Permeability: " << permeability << std::endl;

        // ==================== AUTOMATIC DELTA_P ADJUSTMENT ====================
        // Convert permeability to metric units
        T permeabilityMetric = permeability * dx * dx;
        pcout << "Permeability in metric units: " << permeabilityMetric << " m^2" << std::endl;

        // Calculate suitable deltaP based on target Péclet number: dP = Pe * mu * D / k
        T viscosity = 1e-3;  // Dynamic viscosity (Pa·s)
        T diffusionCoefficient = vec_solute_poreD[0];
        deltaP = Pe * viscosity * diffusionCoefficient / permeabilityMetric;
        pcout << "\n========== AUTO DELTA_P CALCULATOR ==========\n";
        pcout << "Target Pe: " << Pe << std::endl;
        pcout << "Diffusion coefficient: " << diffusionCoefficient << " m^2/s" << std::endl;
        pcout << "Adjusted deltaP: " << deltaP << " Pa" << std::endl;

        // Rerun NS with adjusted deltaP
        NSdomainSetup(nsLattice, createLocalBoundaryCondition3D<T, DESCRIPTOR>(), geometry, deltaP, nsLatticeOmega, pore_dynamics, bounce_back, no_dynamics);

        pcout << "\n========== LBM NS simulation with adjusted deltaP ==========\n";
        ns_convg1.resetValues();
        for (iT = 0; iT < ns_maxiTer_1; ++iT) {
            nsLattice.collideAndStream();
            ns_convg1.takeValue(getStoredAverageEnergy(nsLattice), true);
            if (ns_convg1.hasConverged()) { break; }
        }

        // Recalculate velocities with adjusted deltaP
        PoreMeanU = computeAverage(*computeVelocityNorm(nsLattice, Box3D(1, nx-2, 0, ny-1, 0, nz-1)));
        PoreMaxUx = computeMax(*computeVelocityComponent(nsLattice, Box3D(1, nx-2, 0, ny-1, 0, nz-1), 0));
        DarcyOutletUx = computeAverage(*computeVelocityComponent(nsLattice, Box3D(nx-2, nx-2, 0, ny-1, 0, nz-1), 0));
        DarcyMiddleUx = computeAverage(*computeVelocityComponent(nsLattice, Box3D((nx-1)/2, (nx-1)/2, 0, ny-1, 0, nz-1), 0));
        DarcyInletUx = computeAverage(*computeVelocityComponent(nsLattice, Box3D(1, 1, 0, ny-1, 0, nz-1), 0));

        permeability = computePermeability(nsLattice, nsNu, deltaP, nsLattice.getBoundingBox());
        pcout << "Final Permeability: " << permeability << std::endl;
        pcout << "Final PoreMeanU: " << PoreMeanU << std::endl;
        pcout << "Final PoreMaxUx: " << PoreMaxUx << std::endl;
        pcout << "==============================================\n";
        // ==================== END AUTOMATIC DELTA_P ADJUSTMENT ====================

        if (track_performance == 0) {
            pcout << "Writing velocity VTI... \n";
            writeNsVTI(nsLattice, iT, "nsLatticeFinal_");
            pcout << "Writing checkpoint... \n";
            saveBinaryBlock(nsLattice, str_outputDir + ns_filename + ".chk");
        }

        global::timer("NS").stop();

        T nstime = global::timer("NS").getTime();

        if (ade_maxiTer == 0) {
            pcout << "ade_max_iTer is set to 0. Terminating the simulation.\n";
            return 0;
        }


    /*  =================================== rxn Lattice setup  ===================================  */

    T refNu;
    T refTau;

    if (Pe > thrd) {
        refNu = PoreMeanU * charcs_length / Pe; 
        refTau = refNu * BGK<T>::invCs2 + 0.5;
        if (refTau > 2) {
            pcout << "Reference relaxation time is > 2 (refTau = " << refTau << "). Consider reducing it for numerical accuracy by reducing average flow velocity (e.g. reduce delta_P).\n";
            return -1;
        }
        else if (refTau <= 0.5) {
            pcout << "Reference relaxation time does not satisfy a necessary stability condition for the BGK operator. (tau must be > 0.5, but refTau = " << refTau << ").\n";
            pcout << "Consider increasing average flow velocity (e.g. increase delta_P).\n";
            return -1;
        }
    }
    else {
        refTau = tau;
        refNu = BGK<T>::cs2 * (refTau - 0.5);
    }

    T refOmega = 1/refTau;
    T ade_dt = refNu * dx * dx / vec_solute_poreD[0];
    pcout << "dt = " << ade_dt << ".\n";

    // Calculate mean velocity in metric units and suggested iterations
    T meanUMetric = PoreMeanU * dx / ade_dt;
    T domainLength = nx * dx;
    T flowDuration = domainLength / meanUMetric;
    plint suggestedMaxIter = static_cast<plint>(flowDuration / ade_dt);
    pcout << "\n========== SIMULATION PARAMETERS ==========\n";
    pcout << "Mean velocity (metric): " << meanUMetric << " m/s" << std::endl;
    pcout << "Domain length: " << domainLength << " m" << std::endl;
    pcout << "Flow duration: " << flowDuration << " s" << std::endl;
    pcout << "Suggested ade_maxiTer: " << suggestedMaxIter << std::endl;
    pcout << "Actual Pe: " << meanUMetric * charcs_length / vec_solute_poreD[0] << std::endl;
    pcout << "============================================\n";

    std::vector<T> substrNUinPore(num_of_substrates), substrTAUinPore(num_of_substrates), substrOMEGAinPore(num_of_substrates);
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        if (iS == 0) {
            substrNUinPore[iS]=refNu; substrTAUinPore[iS]=refTau; substrOMEGAinPore[iS]=refOmega;
        }
        else {
            substrNUinPore[iS]   = substrNUinPore[0]*vec_solute_poreD[iS]/vec_solute_poreD[0];
            substrTAUinPore[iS]  = substrNUinPore[iS]*BGK<T>::invCs2+0.5;
            substrOMEGAinPore[iS]= 1/substrTAUinPore[iS];
        }
       
    }

    pcout << "\nInitializing the reaction lattices... \n";
    pcout << "substrTAUinPore = ";
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        pcout << substrTAUinPore[iS] << " ";
    }
    pcout << std::endl;

    if (Pe > thrd) {
        pcout << "Peclet Number (meanU) = " << PoreMeanU * charcs_length / refNu << ", Grid Peclet Number (maxU) = " << PoreMaxUx / refNu << std::endl;
    }
    pcout << "ade_dt = " << ade_dt <<  " s/ts" << std::endl;

    MultiBlockLattice3D<T, BGK> substrLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T, BGK>(refOmega));

    pcout << "=== substrLattice ===" << std::endl;
    pcout << "Dimensions: " << nx << "x" << ny << "x" << nz << std::endl;
    pcout << "Reference relaxation frequency (refOmega): " << refOmega << std::endl;

    pcout << "Before Initializing vec_substr_lattices" << std::endl;
    std::vector< MultiBlockLattice3D<T, BGK> > vec_substr_lattices(num_of_substrates, substrLattice);
    pcout << "After Initializing vec_substr_lattices with " << num_of_substrates << " substrates." << std::endl;

    std::vector< MultiBlockLattice3D<T, BGK> > dC(num_of_substrates, substrLattice);
    std::vector< MultiBlockLattice3D<T, BGK> > dC0(num_of_substrates, substrLattice);

    pcout << "=== Setting Up Lattices for Concentrations and Concentration Gradients ===" << std::endl;
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        pcout << "Setting up lattice for substrate: " << iS << std::endl;

        soluteDomainSetup(vec_substr_lattices[iS], createLocalAdvectionDiffusionBoundaryCondition3D<T, BGK>(), geometry, substrOMEGAinPore[iS], pore_dynamics, bounce_back, no_dynamics,
            vec_c0[iS], vec_left_btype[iS], vec_right_btype[iS], vec_left_bcondition[iS], vec_right_bcondition[iS]);

        soluteDomainSetup(dC[iS], createLocalAdvectionDiffusionBoundaryCondition3D<T, BGK>(), geometry, substrOMEGAinPore[iS], pore_dynamics, bounce_back, no_dynamics,
            0., vec_left_btype[iS], vec_right_btype[iS], vec_left_bcondition[iS], vec_right_bcondition[iS]);
    }

    pcout << "Before Copying dC to dC0" << std::endl;
    dC0 = dC;
    pcout << "After Copying dC to dC0" << std::endl;

    plint tmpIT0 = 0, tmpIT1 = 0;

    std::vector<plint> loctrack;

    pcout << "Before Initializing  MultiBlockLattice object to store the material numbers for each cell in the lattice." << std::endl;
    MultiBlockLattice3D<T, BGK> maskLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T, BGK>(0.));
    pcout << "after Initializing  MultiBlockLattice object to store the material numbers for each cell in the lattice." << std::endl;

    pcout << "Before Initializing  MultiBlockLattice object to store the distance from the solid surface for each cell in the lattice." << std::endl;
    MultiBlockLattice3D<T, BGK> distLattice(nx, ny, nz, new AdvectionDiffusionBGKdynamics<T, BGK>(0.));
    pcout << "after Initializing  MultiBlockLattice object to store the distance from the solid surface for each cell in the lattice." << std::endl;

    pcout << "Before  defineMaskLatticeDynamics" << std::endl;
    defineMaskLatticeDynamics( substrLattice, maskLattice, thrd );
    pcout << "after  defineMaskLatticeDynamics" << std::endl;

    pcout << "before  applyProcessingFunctional(new CopyGeometryScalar2maskLattice3D" << std::endl;
    applyProcessingFunctional(new CopyGeometryScalar2maskLattice3D<T,BGK,int>(pore_dynamics), maskLattice.getBoundingBox(), maskLattice, geometry);
    pcout << "after  applyProcessingFunctional(new CopyGeometryScalar2maskLattice3D" << std::endl;

    pcout << "Before  applyProcessingFunctional(new CopyGeometryScalar2distLattice3D" << std::endl;
    applyProcessingFunctional(new CopyGeometryScalar2distLattice3D<T, BGK, int>(), distLattice.getBoundingBox(), distLattice, distanceDomain);
    pcout << "after  applyProcessingFunctional(new CopyGeometryScalar2distLattice3D" << std::endl;

    pcout << "Before  pointers set" << std::endl;
    std::vector< MultiBlockLattice3D<T, BGK>* > substrate_lattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        substrate_lattices.push_back(&vec_substr_lattices[iS]);
    }
    substrate_lattices.push_back(&maskLattice);

    std::vector< MultiBlockLattice3D<T, BGK>* > ptr_kns_lattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        ptr_kns_lattices.push_back(&vec_substr_lattices[iS]);
    }
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        ptr_kns_lattices.push_back(&dC[iS]);
    }
    ptr_kns_lattices.push_back(&maskLattice);

    std::vector< MultiBlockLattice3D<T, BGK>* > ptr_update_rxnLattices;
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        ptr_update_rxnLattices.push_back(&vec_substr_lattices[iS]);
    }
    for (plint iS = 0; iS < num_of_substrates; ++iS) {
        ptr_update_rxnLattices.push_back(&dC[iS]);
    }
    ptr_update_rxnLattices.push_back(&maskLattice);
    pcout << "after  pointers set" << std::endl;

    pcout << "solver_type = lbm ";
    for (size_t iT = 0; iT < solver_type.size(); ++iT) {
        if (solver_type[iT] == 1) { pcout << "lbm "; }
    }
    pcout << std::endl;

    for (plint iT2 = 0; iT2 <ns_maxiTer_1; ++iT2) {
        nsLattice.collideAndStream();
        ns_convg1.takeValue(getStoredAverageEnergy(nsLattice),false);
        if (ns_convg1.hasConverged()) { break; }
    }

    if (read_NS_file == 0 || (read_NS_file == 1 && ns_rerun_iT0 > 0)) {
        pcout << std::endl << "Flow calculation finished at iteration = " << iT << std::endl;
        if (track_performance == 0) {
            pcout << "Writing velocity VTI... \n";
            writeNsVTI(nsLattice, ns_maxiTer_1, "nsLatticeFinal1_");
            pcout << "Writing checkpoint... \n";
            saveBinaryBlock(nsLattice, str_outputDir + ns_filename + ".chk");
        }
    }

    if (track_performance == 1) { nstime += global::timer("NS").getTime(); global::timer("NS").stop(); }
    
    if (Pe > thrd) {
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            latticeToPassiveAdvDiff(nsLattice, vec_substr_lattices[iS], vec_substr_lattices[iS].getBoundingBox());  
        }

        pcout << "Stabilizing the ADE lattices after coupling the NS lattice...\n";

        for (plint iT = 0; iT < 10000; ++iT) {
            for (plint iS = 0; iS < num_of_substrates; ++iS) {
                vec_substr_lattices[iS].collideAndStream();  
            }
        }

        for (plint iS = 0; iS < num_of_substrates; ++iS) { 
            applyProcessingFunctional(new stabilizeADElattice<T, BGK, int>(vec_c0[iS], pore_dynamics), vec_substr_lattices[iS].getBoundingBox(), vec_substr_lattices[iS], geometry);  
        }
    }


    /*  ================================= rxn Lattice main loop  =================================  */

    iT = 0;
    
    if (read_ADE_file == 1) {
        if (ade_rerun_iT0 > 0) {
            pcout << "read binary ADE files" << std::endl;
            for (plint iS = 0; iS < num_of_substrates; ++iS) { 
                loadBinaryBlock(vec_substr_lattices[iS], str_outputDir + ade_filename + "_" + std::to_string(iS));
            }
            tmpIT0 = 0; tmpIT1 = 0;
            iT = ade_rerun_iT0;
            pcout << "ADE binary files successfully loaded" << std::endl;
        }
        else {
            pcout << "WARNING: number of input files should be equal to the number of substrates " << std::endl;
            return -1;
        }
    }

    T catime = 0, adetime = 0, knstime = 0, cnstime = 0;

    /*  ================================= Reaction System  =================================  */

    /*  =====================================================================================  */

    pcout << "\n===================== LBM ADE simulation begins =====================\n\n";

    global::timer("ade").restart();

    util::ValueTracer<T> ns_convg2(1.0, 1000.0, ns_converge_iT2);
    
    bool ns_saturate = 0, percolationFlag = 0;

    for (; iT < ade_maxiTer; ++iT) {

        // ========================= save VTI files ========================= //

        if (ade_VTI_iTer > 0 && iT % ade_VTI_iTer == 0) {
            pcout << "Iteration = " << iT << "; current_simulation_time = " << iT * ade_dt << " seconds" << std::endl;
            if (track_performance == 0) {
                for (plint iS = 0; iS < num_of_substrates; ++iS) {
                    writeAdvVTI(vec_substr_lattices[iS], iT, ade_filename + std::to_string(iS) + "_");
                }

                if (Pe > thrd) { writeNsVTI(nsLattice, iT, "nsLattice_"); }
                writeAdvVTI(maskLattice, iT, mask_filename + "_");
                pcout << "Writing ADE VTI files... \n";
            }
            adetime += global::timer("ade").getTime();
            pcout << "(Wall-clock) Time elapsed: " << global::timer("ade").getTime() << " seconds." << std::endl;
            global::timer("ade").restart();
        }

        // ===================== save checkpoint files ====================== //

        if (ade_CHK_iTer > 0 && iT % ade_CHK_iTer == 0 && iT > 0 && track_performance == 0) {

            for (plint iS = 0; iS < num_of_substrates; ++iS) {
                if (vec_fixC.size() > 0) { 
                    if (vec_fixC[iS] == 0) { 
                        saveBinaryBlock(vec_substr_lattices[iS], str_outputDir + ade_filename + std::to_string(iS) + "_" + std::to_string(iT) + ".chk"); 
                    } 
                }
                else { 
                    saveBinaryBlock(vec_substr_lattices[iS], str_outputDir + ade_filename + std::to_string(iS) + "_" + std::to_string(iT) + ".chk"); 
                }
            }

            pcout << "Writing checkpoint files... \n";
        }

        // =================== solute lattice collision ==================== //

        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            vec_substr_lattices[iS].collide();
        }
        
        // ============================ run reaction term ============================ //

        dC = dC0;

        applyProcessingFunctional(new run_kinetics<T, BGK>(nx, num_of_substrates, ade_dt, vec_Kc_kns, no_dynamics, bounce_back),
        vec_substr_lattices[0].getBoundingBox(), ptr_kns_lattices);

        // ======================= update flow field and lattice dynamics =======================
       
        if (track_performance == 1) { global::timer("cns").restart(); }
        for (plint iS = 0; iS < num_of_substrates; ++iS) {
            vec_substr_lattices[iS].stream();
        }
        if (track_performance == 1) { nstime += global::timer("cns").getTime(); global::timer("cns").stop(); }
        if (percolationFlag == 1) { break; }
        pcout << "End of simulation at iteration " << iT << std::endl << std::endl;
       
    }


    /*  ================================= Finalize the simulation  =================================  */

        pcout << "Before Applying  Finalize the simulation" << std::endl;
        T TET = global::timer("total").getTime(); global::timer("total").stop();
        pcout << "Total elapsed time: " << TET << " seconds, " << TET / 60 << " minutes, and " << TET / 3600 << " hours." << std::endl;

        if (track_performance == 1) {
            adetime += global::timer("ade").getTime(); global::timer("ade").stop();
            pcout << "Total time consumed by NS: " << nstime << " seconds, " << nstime / 60 << " minutes, and " << nstime / 3600 << " hours." << std::endl;
            pcout << "Total time consumed by ADE: " << adetime << " seconds, " << adetime / 60 << " minutes, and " << adetime / 3600 << " hours." << std::endl;
            pcout << "Total time consumed by C&S: " << cnstime << " seconds, " << cnstime / 60 << " minutes, and " << cnstime / 3600 << " hours." << std::endl;
            if (kns_count > 0) { pcout << "Total time consumed by KNS: " << knstime << " seconds, " << knstime / 60 << " minutes, and " << knstime / 3600 << " hours." << std::endl; }
        }

        else {
            pcout << "Writing VTI and CHK files ..." << std::endl << std::endl;
            for (plint iS = 0; iS < num_of_substrates; ++iS) {
                writeAdvVTI(vec_substr_lattices[iS], iT, ade_filename + std::to_string(iS) + "_");
                saveBinaryBlock(vec_substr_lattices[iS], str_outputDir + ade_filename + std::to_string(iS) + "_" + std::to_string(iT) + ".chk");
            }
            writeAdvVTI(maskLattice, iT, mask_filename + "_");
            saveBinaryBlock(maskLattice, str_outputDir + mask_filename + "_" + std::to_string(iT) + ".chk");
            if (Pe > thrd) {
                writeNsVTI(nsLattice, iT, "nsLattice_");
                saveBinaryBlock(nsLattice, str_outputDir + ns_filename + ".chk");
            }
        }

        pcout << "after Applying  Finalize the simulation" << std::endl;

    pcout << "\nSimulation Finished!" << std::endl << std::endl;

    return 0;

}
