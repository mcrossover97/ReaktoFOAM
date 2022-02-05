// Reaction Rate Functions

auto theMicroclineRateFunction(const Reaktoro::Reaction& reaction, const Reaktoro::ChemicalSystem& system) -> Reaktoro::ReactionRateFunction
{
    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();
	const unsigned MicroclineIndex = system.indexSpecies("Microcline");
	const unsigned protonIndex = system.indexSpecies("H+");
	const unsigned HFIndex = system.indexSpecies("HF@");

    // Auxiliary variables
    Reaktoro::ChemicalScalar mechanism, g;

    // Define the mineral mechanism function
    Reaktoro::ReactionRateFunction fn = [=](const Reaktoro::ChemicalProperties& properties) mutable
    {
		double specificSurfaceArea = 3100; // cm2/g
		double p = 1;
		double q = 1;
		double ProtonConc = properties.composition()[protonIndex].val / properties.fluidVolume().val / 1000;
		double HFConc = properties.composition()[HFIndex].val / properties.fluidVolume().val / 1000;
		double alpha = 1.2;
		double beta = 0.4;
		double Ef = 0.127 * std::exp(-4680/temp);
		double K = 5.66E-2 * std::exp(956/temp);
		double kappa = Ef * (1 + K * std::pow(ProtonConc, beta)) * std::pow(HFConc, alpha) * 1000;
		double molar_mass = system.species(MicroclineIndex).molarMass();

        // The result of this function evaluation
        Reaktoro::ChemicalScalar res(num_species);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega);

        // Calculate the p and q powers of the saturation index Omega
        const auto pOmega = pow(Omega, p);
        const auto qOmega = pow(1 - pOmega, q);

        // Calculate the function mechanism
        mechanism = kappa * qOmega;
		
        // Calculate the function g
		Reaktoro::ChemicalScalar f (num_species);
		
		f = 0.0;
		
		const auto n = properties.composition();
		auto nm = n[MicroclineIndex];
		
		nm.val = std::max(nm.val, 0.0);
		
		f += mechanism;
		f *= molar_mass * specificSurfaceArea * 1E-1;

        // Calculate the resulting mechanism function
        res = nm * f;

        return res;
    };
	return fn;
}
	
auto theAlbiteRateFunction(const Reaktoro::Reaction& reaction, const Reaktoro::ChemicalSystem& system) -> Reaktoro::ReactionRateFunction
{
    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();
	const unsigned AlbiteIndex = system.indexSpecies("Albite");
	const unsigned protonIndex = system.indexSpecies("H+");
	const unsigned HFIndex = system.indexSpecies("HF@");

    // Auxiliary variables
    Reaktoro::ChemicalScalar mechanism, g;

    // Define the mineral mechanism function
    Reaktoro::ReactionRateFunction fn = [=](const Reaktoro::ChemicalProperties& properties) mutable
    {
		double specificSurfaceArea = 200; // cm2/g
		double p = 1;
		double q = 1;
		double ProtonConc = properties.composition()[protonIndex].val / properties.fluidVolume().val / 1000;
		double HFConc = properties.composition()[HFIndex].val / properties.fluidVolume().val / 1000;
		double alpha = 1;
		double beta = 1;
		double Ef = 9.5E-3 * std::exp(-3930/temp);
		double K = 6.24E-2 * std::exp(554/temp);
		double kappa = Ef * (1 + (K * std::pow(ProtonConc, beta))) * std::pow(HFConc, alpha) * 1000;
		double molar_mass = system.species(AlbiteIndex).molarMass();

        // The result of this function evaluation
        Reaktoro::ChemicalScalar res(num_species);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega);

        // Calculate the p and q powers of the saturation index Omega
        const auto pOmega = pow(Omega, p);
        const auto qOmega = pow(1 - pOmega, q);

        // Calculate the function mechanism
        mechanism = kappa * qOmega;
		
        // Calculate the function g
		Reaktoro::ChemicalScalar f (num_species);
		
		f = 0.0;
		
		const auto n = properties.composition();
		auto nm = n[AlbiteIndex];
		
		nm.val = std::max(nm.val, 0.0);
		
		f += mechanism;
		f *= molar_mass * specificSurfaceArea * 1E-1;

        // Calculate the resulting mechanism function
        res = nm * f;

        return res;
    };
	return fn;
}

auto theQuartzRateFunction(const Reaktoro::Reaction& reaction, const Reaktoro::ChemicalSystem& system) -> Reaktoro::ReactionRateFunction
{
    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();
	const unsigned QuartzIndex = system.indexSpecies("Quartz");
	const unsigned HFIndex = system.indexSpecies("HF@");

    // Auxiliary variables
    Reaktoro::ChemicalScalar mechanism, g;

    // Define the mineral mechanism function
    Reaktoro::ReactionRateFunction fn = [=](const Reaktoro::ChemicalProperties& properties) mutable
    {
		double specificSurfaceArea = 200; // cm2/g
		double p = 1;
		double q = 1;
		double HFConc = properties.composition()[HFIndex].val / properties.fluidVolume().val / 1000;
		double Ef = 1.39E-7 * std::exp(-1150/temp);
		double kappa = Ef * HFConc * 1000;
		double molar_mass = system.species(QuartzIndex).molarMass();

        // The result of this function evaluation
        Reaktoro::ChemicalScalar res(num_species);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega);

        // Calculate the p and q powers of the saturation index Omega
        const auto pOmega = pow(Omega, p);
        const auto qOmega = pow(1 - pOmega, q);

        // Calculate the function mechanism
        mechanism = kappa * qOmega;
		
        // Calculate the function g
		Reaktoro::ChemicalScalar f (num_species);
		
		f = 0.0;
		
		const auto n = properties.composition();
		auto nm = n[QuartzIndex];
		
		nm.val = std::max(nm.val, 0.0);
		
		f += mechanism;
		f *= molar_mass * specificSurfaceArea * 1E-1;

        // Calculate the resulting mechanism function
        res = nm * f;

        return res;
    };
	return fn;
}

auto theKaoliniteRateFunction(const Reaktoro::Reaction& reaction, const Reaktoro::ChemicalSystem& system) -> Reaktoro::ReactionRateFunction
{
    // The number of chemical species in the system
    const unsigned num_species = system.numSpecies();
	const unsigned KaoliniteIndex = system.indexSpecies("Kaolinite");
	const unsigned HFIndex = system.indexSpecies("HF@");

    // Auxiliary variables
    Reaktoro::ChemicalScalar mechanism, g;

    // Define the mineral mechanism function
    Reaktoro::ReactionRateFunction fn = [=](const Reaktoro::ChemicalProperties& properties) mutable
    {
		double specificSurfaceArea = 20000; // cm2/g
		double p = 1;
		double q = 1;
		double HFConc = properties.composition()[HFIndex].val / properties.fluidVolume().val / 1000;
		double alpha = 1;
		double Ef = 0.86 * std::exp(-6800/temp);
		double kappa = Ef * std::pow(HFConc, alpha) * 1000;
		double molar_mass = system.species(KaoliniteIndex).molarMass();

        // The result of this function evaluation
        Reaktoro::ChemicalScalar res(num_species);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        // Calculate the saturation index
        const auto Omega = exp(lnOmega);

        // Calculate the p and q powers of the saturation index Omega
        const auto pOmega = pow(Omega, p);
        const auto qOmega = pow(1 - pOmega, q);

        // Calculate the function mechanism
        mechanism = kappa * qOmega;
		
        // Calculate the function g
		Reaktoro::ChemicalScalar f (num_species);
		
		f = 0.0;
		
		const auto n = properties.composition();
		auto nm = n[KaoliniteIndex];
		
		nm.val = std::max(nm.val, 0.0);
		
		f += mechanism;
		f *= molar_mass * specificSurfaceArea * 1E-1;

        // Calculate the resulting mechanism function
        res = nm * f;

        return res;
    };
	return fn;
}