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
    PLEAFOAM

Description
   A Geochemical Reactive Transport Solver.
	
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <Reaktoro/Reaktoro.hpp>
#include <ThermoFun/ThermoFun.h>
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double temp = 298.15;
#include <ReactionRateFunctions.h>

int main(int argc, char *argv[])
{	
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
	#include "createTimeControls.H"
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	
	// Defining some initial properties
	int P0 = 101325;
	double T = temp.value();
	double refBulkVol = 1;
	double inletCellNo = 0;
	double columnSum;
	double averagePermRatio;
	
	// Defining patches
	label patchi = mesh.boundaryMesh().findPatchID("inlet");
	const polyPatch& cPatch = mesh.boundaryMesh()[patchi];
	label patchii = mesh.boundaryMesh().findPatchID("outlet");	
	
	Info<< "\ninitialization starts...\n" << endl;
	
	ThermoFun::Database database("databases/mines16-thermofun.json");
	
	std::string theElements = "H O F Cl Si Al K Na";
	
	Reaktoro::ChemicalEditor editor(database);
	editor.setTemperatures({T}, "celsius");
	editor.addAqueousPhaseWithElements(theElements);
	//	.setChemicalModelDebyeHuckel();
		
	editor.addMineralPhase("Quartz");
	editor.addMineralPhase("Albite");
	editor.addMineralPhase("Microcline");
	editor.addMineralPhase("Kaolinite");
	
	editor.addMineralPhase("AluminumFluoride");
	editor.addMineralPhase("Analcime");
	editor.addMineralPhase("Boehmite");
	editor.addMineralPhase("Chabazite-(Na)");
	editor.addMineralPhase("Chalcedony");
	editor.addMineralPhase("Coesite");
	editor.addMineralPhase("Corundum");
	editor.addMineralPhase("Dickite");
	editor.addMineralPhase("Gibbsite");
	editor.addMineralPhase("Halite");
	editor.addMineralPhase("Heulandite-(Na)");
	editor.addMineralPhase("Hi-Albite");
	editor.addMineralPhase("Jadeite");
	editor.addMineralPhase("K-Beidellite");
	editor.addMineralPhase("Kalsilite");
	editor.addMineralPhase("Kyanite");
	editor.addMineralPhase("Leucite");
	editor.addMineralPhase("Mordenite-(Na)");
	editor.addMineralPhase("Na-Beidellite");
	editor.addMineralPhase("Natrolite");
	editor.addMineralPhase("Nepheline");
	editor.addMineralPhase("PotassiumFluorosilicate");
	editor.addMineralPhase("Pyrophyllite");
	editor.addMineralPhase("Sanidine");
	editor.addMineralPhase("Sillimanite");
	editor.addMineralPhase("SodiumFluorosilicate");
	editor.addMineralPhase("SodiumHexafluoroaluminate");
	editor.addMineralPhase("SodiumOxide");
	editor.addMineralPhase("Sylvite");
	editor.addMineralPhase("Topaz-F");
	editor.addMineralPhase("Topaz-OH");
	editor.addMineralPhase("Tridymite");
	
	//editor.addMineralPhase("Muscovite");
	//editor.addMineralPhase("Paragonite");

	Reaktoro::ChemicalSystem system(editor);
	
	int c = 0;
	for(auto species : system.species()){
		std::cout<<species.name()<<std::endl;
		c++;
	}
	
	
	#include <reactions.h>
	Reaktoro::ReactionSystem reactionSystem(system, reactions);
	Reaktoro::Partition partition(system);
	const std::vector<std::string> kinetikSpecies{"Quartz","Microcline","Albite","Kaolinite"};
	partition.setKineticSpecies(kinetikSpecies);
	
	Reaktoro::KineticPath path(reactionSystem);
    path.setPartition(partition);
	
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
	
	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\ninitializing inlet...\n" << endl;
	
	Reaktoro::EquilibriumProblem problem_bc(system);
	problem_bc.setTemperature(T,"celsius");
	problem_bc.setPressure(P0);
	problem_bc.add("H2O", 85, "g");
	problem_bc.add("HCl", 12, "g");
	problem_bc.add("HF", 3, "g");
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
		state_ic.scalePhaseVolume("Quartz", (1-eps[celli])*refBulkVol*0.75, "m3");
		state_ic.scalePhaseVolume("Microcline", (1-eps[celli])*refBulkVol*0.05, "m3");
		state_ic.scalePhaseVolume("Albite", (1-eps[celli])*refBulkVol*0.05, "m3");
		state_ic.scalePhaseVolume("Kaolinite", (1-eps[celli])*refBulkVol*0.15, "m3");
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
		
		Info<< "\nEquilibrate" << nl << endl;
		forAll(states,celli){	
		
			Reaktoro::EquilibriumProblem problemRun(system);
			problemRun.setPartition(partition);
			problemRun.setTemperature(T, "celsius");
			problemRun.setPressure(p[celli]);
			
			i = 0;
			for(auto elements : system.elements()){
				problemRun.add(elements.name(), FluidElements[i][celli] , "mol");
				i++;
			}	
			
			solver.solve(states[celli],problemRun);
			i = 0;
			j = 0;
			
			for(auto species : system.species()){
				if (i<numberOfAqueousSpecies){
					i++;
				} else {
					states[celli].setSpeciesAmount(species.name(),SolidSpecies[j][celli], "mol");
					j++;
				}
			}
			
			path.solve(states[celli], 0, runTime.deltaTValue() , "second");
			
			for(auto species : system.species()){
				if(states[celli].speciesAmount(species.name()) < 1e-20){
					states[celli].setSpeciesAmount(species.name(), 1e-20, "mol") ;
				}
			}
			
			properties[celli] = states[celli].properties();
			eps[celli] = 1.0-properties[celli].solidVolume().val;
			k[celli] = k0.value()*pow(eps[celli]/eps0.value(),3);
				
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
		
		//Info<< "Calculating Velocities ..." << nl << endl;
		//U = -k/mu*fvc::grad(p);
		//
		//Info<< "Correcting BCs ..." << nl << endl;
		//U.correctBoundaryConditions();
		
		Info<< "Calculating Fluxes ..." << nl << endl;
		phi = - pEqn.flux();
		
		Info<< "updating bc pressures" << nl << endl;
		forAll(cPatch,facei){
			Reaktoro::EquilibriumProblem problem_bc(system);
			problem_bc.setTemperature(T,"celsius");
			problem_bc.setPressure(p.boundaryField()[patchi][facei]);
			problem_bc.add("H2O", 85, "g");
			problem_bc.add("HCl", 12, "g");
			problem_bc.add("HF", 3, "g");
			Reaktoro::ChemicalState state_bc = Reaktoro::equilibrate(problem_bc);
			state_bc.scaleVolume(refBulkVol, "m3");
			inletStates.set(facei, new Reaktoro::ChemicalState(state_bc));
		}
			
		Info<< "Solving Fluid Elements Mass Balance Equation ..." << nl << endl;
		for (int i = 0; i<numberOfTotalElements; i++){
			fvScalarMatrix transportEqn
			(
				fvm::ddt(1.0,FluidElements[i]) + fvm::div(phi/fvc::interpolate(eps),FluidElements[i]) -  fvm::laplacian(De,FluidElements[i])
			);
			transportEqn.solve();	
		}
		
		runTime.write();
		
		//columnSum=0;
		//for(int i=0; i<xNo.value(); i++)
		//{
		//	columnSum = columnSum + 1/k[i];
		//}
		//
		//averagePermRatio=xNo.value()/columnSum/9.869233e-16;
		//
		//ofstream out("AveragePermeability", ios::app);
		//out << averagePermRatio << endl;
		//
		//ofstream tout("Time", ios::app);
		//tout << runTime.timeName() << endl;
		
		forAll(cPatch,facei){
			ofstream out("AveragePermeability", ios::app);
			out << 0.00015257843*5*0.001*0.0508/(p.boundaryFieldRef()[patchi][facei]-101325)/9.869233e-16 << endl;
			
			ofstream tout("Time", ios::app);
			tout << runTime.timeName() << endl;
		}
		
		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< nl << endl;
    }
	
    Info<< "End\n" << endl;
	runTime.writeNow();
    return 0;
}


// ************************************************************************* //
