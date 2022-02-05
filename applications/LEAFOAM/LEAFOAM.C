/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    LEAFOAM

Description
   A Geochemical Reactive Transport Solver using LEA.
	
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <Reaktoro/Reaktoro.hpp>
#include <ThermoFun/ThermoFun.h>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{	
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
	#include "createTimeControls.H"
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	int P0 = 101325;
	double T = temp.value();
	double refBulkVol = 1;
	double inletCellNo = 0;
	double ttime = 0;
	double sum = 0;
	double porro = 3.1;
	
	label patchi = mesh.boundaryMesh().findPatchID("inlet");
	const polyPatch& cPatch = mesh.boundaryMesh()[patchi];
	
	label patchii = mesh.boundaryMesh().findPatchID("outlet");	
	
	Info<< "\ninitialization starts...\n" << endl;
	
	ThermoFun::Database database("databases/mines16-thermofun-original.json");
	
	std::string theElements = "H O Cl Si Al C Ca K Mg Na Fe S";
	Reaktoro::ChemicalEditor editorInit(database);
	editorInit.setTemperatures({T}, "celsius");
	editorInit.addAqueousPhaseWithElements(theElements)
		.setChemicalModelHKF();
	editorInit.addMineralPhaseWithElements(theElements);
	Reaktoro::ChemicalSystem systemInit(editorInit);
	int numberOfAqueousSpeciesInit = systemInit.numSpeciesInPhase (0);
	
	Reaktoro::ChemicalEditor editor(database);
	editor.setTemperatures({T}, "celsius");

	editor.addAqueousPhase({"H2O@","OH-","H+","NaCl@","Na+","Cl-","HCl@","SiO2@","H2@","H2O2@","NaOH@","K+","Mg+2","Al+3","Ca+2","CO3-2","HS-","H2S@","CO2@","CaCl2@","Ca(CO3)@","Al(OH)2+","Al(OH)4-","CaCl+","CaOH+","KOH@","KCl@","HCO3-","Fe+2","Fe+3","FeCl+","FeCl2@","Na(CO3)-","Mg(CO3)@","Mg(HCO3)+","MgCl+","MgOH+","AlH3SiO4+2","AlOH+2","Ca(HCO3)+","Ca(HSiO3)+","KAlO2@","Na(HCO3)@","NaAl(OH)4@","NaHSiO3@","Fe(CO3)@","Fe(HCO3)+","FeO2H-","FeO@","FeOH+","Mg(HSiO3)+","S3-2","S4-2","S5-2","Al(OH)3@","Ca(SO4)@","ClO-","ClO4-","FeCl+2","FeCl2+","FeCl3@","FeO+","FeO2-","FeO2H@","FeOH+2","H2S2O3@","H2S2O4@","HClO@","HO2-","HS2O3-","K(SO4)-","MgSO4@","Na(SO4)-","O2@","S2O8-2","S3O6-2","S4O6-2","S5O6-2","S2O3-2","S2O4-2","S2O5-2","S2O6-2","SO4-2","SO3-2","SO2@"})
		.setChemicalModelDebyeHuckel();
	
	//,"CH4@",
	
	int c = 0;
	for(auto species : systemInit.species()){
		std::cout<<species.name()<<std::endl;
		c++;
	}
	
	editor.addMineralPhase("Albite");
	editor.addMineralPhase("Ankerite");
	editor.addMineralPhase("Calcite");
	editor.addMineralPhase("Chlorite-Mg");
	editor.addMineralPhase("Dolomite");
	editor.addMineralPhase("Microcline");
	editor.addMineralPhase("Muscovite");
	editor.addMineralPhase("Pyrite");
	editor.addMineralPhase("Quartz");
	editor.addMineralPhase("Magnetite");
	editor.addMineralPhase("Kaolinite");

	//editor.addMineralPhase("Aegirine");
	//editor.addMineralPhase("Analcime");
	//editor.addMineralPhase("Anhydrite");
	//editor.addMineralPhase("Annite");
	//editor.addMineralPhase("Anorthite");
	//editor.addMineralPhase("Ca-Beidellite");
	//editor.addMineralPhase("Ca-Montmorill");
	//editor.addMineralPhase("Ca-Nontronite");
	//editor.addMineralPhase("Chabazite-(Ca)");
	//editor.addMineralPhase("Chabazite-(Na)");
	//editor.addMineralPhase("Chalcedony");
	//editor.addMineralPhase("Chrysotile");
	//editor.addMineralPhase("Clinochlore");
	//editor.addMineralPhase("Clinozoisite");
	//editor.addMineralPhase("Dickite");
	//editor.addMineralPhase("Diopside");
	//editor.addMineralPhase("Enstatite");
	//editor.addMineralPhase("Epidote");
	//editor.addMineralPhase("Fayalite");
	//editor.addMineralPhase("Fe(OH)3");
	//editor.addMineralPhase("Fe-Epidote");
	//editor.addMineralPhase("Fe-Sudoite");
	//editor.addMineralPhase("Fe-Talc");
	//editor.addMineralPhase("Ferroactinolite");
	//editor.addMineralPhase("Forsterite");
	//editor.addMineralPhase("Goethite");
	//editor.addMineralPhase("Grossular");
	//editor.addMineralPhase("Hematite");
	//editor.addMineralPhase("Heulandite-(Ca)");
	//editor.addMineralPhase("Heulandite-(Na)");
	//editor.addMineralPhase("Hi-Albite");
	//editor.addMineralPhase("K-Beidellite");
	//editor.addMineralPhase("K-Montmorill");
	//editor.addMineralPhase("K-Nontronite");
	//editor.addMineralPhase("Laumontite");
	//editor.addMineralPhase("Leonhardite");
	//editor.addMineralPhase("Magnesite");
	//editor.addMineralPhase("Mesolite");
	//editor.addMineralPhase("Mg-Beidellite");
	//editor.addMineralPhase("Mg-Montmorill");
	//editor.addMineralPhase("Mg-Nontronite");
	//editor.addMineralPhase("Mordenite-(Na)");
	//editor.addMineralPhase("Na-Beidellite");
	//editor.addMineralPhase("Na-Montmorill");
	//editor.addMineralPhase("Na-Nontronite");
	//editor.addMineralPhase("Natrolite");
	//editor.addMineralPhase("Nepheline");
	//editor.addMineralPhase("Prehnite");
	//editor.addMineralPhase("Pumpellyite");
	//editor.addMineralPhase("Pyrophyllite");
	//editor.addMineralPhase("Pyrrhotite");
	//editor.addMineralPhase("Sanidine");
	//editor.addMineralPhase("Scolecite");
	//editor.addMineralPhase("Siderite");
	//editor.addMineralPhase("Stilbite-(Ca)");
	//editor.addMineralPhase("Thompsonite");
	//editor.addMineralPhase("Topaz-OH");
	//editor.addMineralPhase("Wairakite");
	//editor.addMineralPhase("Wollastonite");
	//editor.addMineralPhase("Yugawaralite");
	
	//editor.addMineralPhase("Arfvedsonite");
	//editor.addMineralPhase("Daphnite");
	//editor.addMineralPhase("Fe-Celadonite");
	//editor.addMineralPhase("Fe-Saponite");
	//editor.addMineralPhase("Mg-Celadonite");
	//editor.addMineralPhase("Riebeckite");
	//editor.addMineralPhase("Amesite-Mg");
	//editor.addMineralPhase("Ca-Saponite");
	//editor.addMineralPhase("Diamond");
	//editor.addMineralPhase("Graphite");
	//editor.addMineralPhase("K-Saponite");
	//editor.addMineralPhase("Mg-Saponite");
	//editor.addMineralPhase("Mordenite-(Ca)");
	//editor.addMineralPhase("Na-Saponite");
	//editor.addMineralPhase("Talc");

	Reaktoro::ChemicalSystem system(editor);
	Reaktoro::EquilibriumSolver solver(system);
	
	int numberOfAqueousSpecies = system.numSpeciesInPhase (0);
	int numberOfTotalSpecies = system.numSpecies();
	int numberOfTotalElements = system.numElements();
	
	PtrList<Reaktoro::ChemicalState> states(mesh.nCells());
	forAll(cPatch,facei){
		inletCellNo ++;
	}
	PtrList<Reaktoro::ChemicalState> inletStates(inletCellNo);
	PtrList<Reaktoro::ChemicalProperties> properties(mesh.nCells());
	PtrList<volScalarField> FluidElements(numberOfTotalElements);
	PtrList<volScalarField> SolidSpecies(numberOfTotalSpecies-numberOfAqueousSpecies);
	PtrList<volScalarField> FluidSpecies(numberOfAqueousSpecies);
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\ninitializing inlet...\n" << endl;
	
	Reaktoro::EquilibriumProblem problem_bc(system);
	problem_bc.setTemperature(T,"celsius");
	problem_bc.setPressure(P0);
	problem_bc.add("H2O", 1, "kg");
	problem_bc.add("NaCl", 0.0005, "mol");
	Reaktoro::ChemicalState state_bc = Reaktoro::equilibrate(problem_bc);
	state_bc.scaleVolume(refBulkVol, "m3");
	forAll(cPatch,facei){
		inletStates.set(facei, new Reaktoro::ChemicalState(state_bc));
	}
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	Info<< "\ninitializing internal field!\n" << endl;
	
	forAll(states,celli){
		Reaktoro::EquilibriumProblem problem_ic(system);
		problem_ic.setTemperature(T, "celsius");
		problem_ic.setPressure(P0);
		problem_ic.add("H2O", 1, "kg");
		Reaktoro::ChemicalState state_ic = Reaktoro::equilibrate(problem_ic);
		state_ic.scalePhaseVolume("Aqueous", eps[celli]*refBulkVol, "m3");
		state_ic.scalePhaseVolume("Calcite", (1-eps[celli])*refBulkVol*0.013, "m3");
		state_ic.scalePhaseVolume("Dolomite", (1-eps[celli])*refBulkVol*0.037/2, "m3");
		state_ic.scalePhaseVolume("Ankerite", (1-eps[celli])*refBulkVol*0.037/2, "m3");
		state_ic.scalePhaseVolume("Pyrite", (1-eps[celli])*refBulkVol*0.17, "m3");
		state_ic.scalePhaseVolume("Quartz", (1-eps[celli])*refBulkVol*0.298, "m3");
		state_ic.scalePhaseVolume("Albite", (1-eps[celli])*refBulkVol*0.044, "m3");
		state_ic.scalePhaseVolume("Microcline", (1-eps[celli])*refBulkVol*0.009, "m3");
		state_ic.scalePhaseVolume("Chlorite-Mg", (1-eps[celli])*refBulkVol*0.099, "m3");
		state_ic.scalePhaseVolume("Muscovite", (1-eps[celli])*refBulkVol*0.33, "m3");
		states.set(celli, new Reaktoro::ChemicalState(state_ic));
		Reaktoro::ChemicalProperties property_ic = state_ic.properties();
		properties.set(celli, new Reaktoro::ChemicalProperties(property_ic));
		eps[celli] = 1.0-properties[celli].solidVolume().val;
		std::cout<<eps[celli]<<std::endl;
	}
	
	int i = 0;
	for(auto elements : system.elements()){
		
		FluidElements.set
		(
			i,
			new volScalarField 
			(
				IOobject
				(
					"Mol_Fluid_"+elements.name(),
					runTime.timeName(),
					mesh,
					IOobject::NO_READ,
					IOobject::NO_WRITE
				),
				mesh,
				dimensionedScalar("Mol_Fluid_"+elements.name(),dimensionSet(0,0,0,0,1,0,0), 0.0)
			)
		);
	
		FluidElements[i].boundaryFieldRef().set(patchi, fvPatchField<scalar>::New("fixedValue", mesh.boundary()[patchi], FluidElements[i]));
		FluidElements[i].boundaryFieldRef().set(patchii, fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchii], FluidElements[i]));

		forAll(mesh.cells(),celli){
			FluidElements[i][celli] = states[celli].elementAmountInPhase(elements.name(),"Aqueous");
		}
		
		forAll(cPatch, faceI){
			FluidElements[i].boundaryFieldRef()[patchi][faceI] = inletStates[faceI].elementAmount(elements.name());
		}
		i++;
	}
	
	i = 0;
	int j = 0;
	for(auto species : system.species()){
		if (i<numberOfAqueousSpecies){
			FluidSpecies.set
			(
				i,
				new volScalarField 
				(
					IOobject
					(
						"Mol_F_"+species.name(),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::NO_WRITE
					),
					mesh,
					dimensionedScalar("Mol_F_"+species.name(),dimensionSet(0,0,0,0,1,0,0), 0.0)
				)
			);
			
			forAll(mesh.cells(),celli){
				FluidSpecies[i][celli] = states[celli].speciesAmount(species.name());
			}
			
			i++;
		}
		else {
			SolidSpecies.set
			(
				j,
				new volScalarField 
				(
					IOobject
					(
						"Mol_Solid_"+species.name(),
						runTime.timeName(),
						mesh,
						IOobject::NO_READ,
						IOobject::AUTO_WRITE
					),
					mesh,
					dimensionedScalar("Mol_Solid_"+species.name(),dimensionSet(0,0,0,0,1,0,0), 0.0)
				)
			);
			
			SolidSpecies[j].boundaryFieldRef().set(patchi, fvPatchField<scalar>::New("fixedValue", mesh.boundary()[patchi], SolidSpecies[j]));
			SolidSpecies[j].boundaryFieldRef().set(patchii, fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchii], SolidSpecies[j]));
			
			forAll(mesh.cells(),celli){
				SolidSpecies[j][celli] = states[celli].speciesAmount(species.name());
			}
			
			forAll(cPatch, faceI){
				SolidSpecies[j].boundaryFieldRef()[patchi][faceI] = 0;
			}
			
			j++;
		}
	}
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
    Info<< "\nCalculation starts!\n" << endl;
	
    while (runTime.loop())
    {
		#include "readTimeControls.H"
		#include "CourantNo.H"
		#include "setDeltaT.H"
		
		Info<< "Time = " << runTime.timeName() << nl << endl;
		
		ttime = runTime.value();
		
		Info<< "\nEquilibrate" << nl << endl;
		forAll(states,celli){
			Reaktoro::EquilibriumProblem problemRun(system);
			problemRun.setTemperature(T, "celsius");
			//problemRun.setPressure(p[celli]);
			
			i = 0;
			for(auto elements : system.elements()){
				problemRun.add(elements.name(), FluidElements[i][celli] , "mol");
				i++;
			}	

			i = 0;
			j = 0;
			
			for(auto species : system.species()){
				if (i<numberOfAqueousSpecies){
					i++;
				} else {
					problemRun.add(species.name(), SolidSpecies[j][celli], "mol");
					j++;
				}
			}	
			
			solver.solve(states[celli],problemRun);
			
			properties[celli] = states[celli].properties();
			eps[celli] = 1.0-properties[celli].solidVolume().val;
			if (eps[celli]<0 || eps[celli]>1){
				try {
					Reaktoro::ChemicalState stat = Reaktoro::equilibrate(problemRun);
					states[celli] = stat;
				} catch (...){
					solver.approximate(states[celli],problemRun);
				}
				properties[celli] = states[celli].properties();
				eps[celli] = 1.0-properties[celli].solidVolume().val;
			} else {
			}
			k[celli] = k0.value()*Foam::pow(eps[celli]/eps0.value(),3.10);
			Info<< eps[celli] << endl;
			
			i = 0;  
			for(auto elements : system.elements()){
				FluidElements[i][celli] = states[celli].elementAmountInPhase(elements.name(),"Aqueous");
				i++;
			}
			
			i = 0;
			j = 0;
			
			for(auto species : system.species()){
				if (i<numberOfAqueousSpecies){
					FluidSpecies[i][celli] = states[celli].speciesAmount(species.name());
					i++;
				} else {
					SolidSpecies[j][celli] = states[celli].speciesAmount(species.name());
					j++;
				}
			}
		}
		
		eps.correctBoundaryConditions();
		k.correctBoundaryConditions();
		
		Info<< "Solving Continuity Equation ..." << nl << endl;
		fvScalarMatrix pEqn
		(
			fvm::laplacian(k/mu,p) - fvc::ddt(1.0,eps)
		);
		pEqn.solve();
		
		Info<< "Calculating Fluxes ..." << nl << endl;
		Info<< ttime << nl << endl;
		phi = - pEqn.flux();
		
		Info<< "updating bc cond" << nl << endl;
		forAll(cPatch,facei){
			if (ttime<805){
				Info<< "yes" << nl << endl;
				//Reaktoro::EquilibriumProblem problem_bc(system);
				//problem_bc.setTemperature(T,"celsius");
				//problem_bc.add("H2O", 1, "kg");
				//problem_bc.add("NaCl", 0.0005, "mol");
				//Reaktoro::ChemicalState state_bc = Reaktoro::equilibrate(problem_bc);
				//state_bc.scaleVolume(refBulkVol, "m3");
				//inletStates.set(facei, new Reaktoro::ChemicalState(state_bc));
				Reaktoro::EquilibriumInverseProblem problem_bcc(system);
				problem_bcc.setTemperature(T,"celsius");
				problem_bcc.add("H2O", 1, "kg");
				problem_bcc.add("NaCl", 0.01, "kg");
				problem_bcc.pH(2.0, "HCl", "NaCl");
				Reaktoro::ChemicalState state_bcc = Reaktoro::equilibrate(problem_bcc);
				state_bcc.scaleVolume(refBulkVol, "m3");
				inletStates.set(facei, new Reaktoro::ChemicalState(state_bcc));
			} else {
				Info<< "no" << nl << endl;
				Reaktoro::EquilibriumInverseProblem problem_bcc(system);
				problem_bcc.setTemperature(T,"celsius");
				problem_bcc.add("H2O", 1, "kg");
				problem_bcc.add("NaCl", 0.01, "kg");
				problem_bcc.pH(2.0, "HCl", "NaCl");
				Reaktoro::ChemicalState state_bcc = Reaktoro::equilibrate(problem_bcc);
				state_bcc.scaleVolume(refBulkVol, "m3");
				inletStates.set(facei, new Reaktoro::ChemicalState(state_bcc));
			}
			//problem_bc.setPressure(p.boundaryField()[patchi][facei]);
		}
		
		i = 0;
		for(auto elements : system.elements()){
			forAll(cPatch, faceI){
				FluidElements[i].boundaryFieldRef()[patchi][faceI] = inletStates[faceI].elementAmount(elements.name());
			}
			i++;
		}
			
		Info<< "Solving Fluid Elements Mass Balance Equation ..." << nl << endl;
		for (int i = 0; i<numberOfTotalElements; i++){
			fvScalarMatrix transportEqn
			(
				fvm::ddt(1.0,FluidElements[i]) + fvm::div(phi/fvc::interpolate(eps),FluidElements[i]) -  fvm::laplacian(De,FluidElements[i])
			);
			transportEqn.solve();	
		}

		sum = 0;
		forAll(eps,cellII){
			sum = sum + 1/eps[cellII];
		}
		porro = 600/sum;
		Info<< porro << nl << endl;
		
		ofstream out("AveragePermeability", ios::app);
		out << 0.00015257843*0.001*0.0508/(p[0]-101325)/9.869233000000002e-19 << endl;
			
		ofstream tout("Time", ios::app);
		tout << runTime.timeName() << endl;
		
		ofstream dout("Porosity", ios::app);
		dout << porro << endl;

		runTime.write();
		
		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;
    }
	
    Info<< "End\n" << endl;
	runTime.writeNow();
    return 0;
}


// ************************************************************************* //
