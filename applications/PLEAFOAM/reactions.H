// Reaction Specifications

const Reaktoro::ReactionEquation equationMicrocline("Microcline + 2*H2O@ = K+ + Al(OH)4- + 3*SiO2@");
Reaktoro::Reaction reactionMicrocline(equationMicrocline, system);
reactionMicrocline.setName("reactionMicrocline");
Reaktoro::ReactionRateFunction ReactionMicroclineRateFunction = theMicroclineRateFunction(reactionMicrocline,system);
reactionMicrocline.setRate(ReactionMicroclineRateFunction);

const Reaktoro::ReactionEquation equationAlbite("Albite + 2*H2O@ = Na+ + Al(OH)4- + 3*SiO2@");
Reaktoro::Reaction reactionAlbite(equationAlbite, system);
reactionAlbite.setName("reactionAlbite");
Reaktoro::ReactionRateFunction ReactionAlbiteRateFunction = theAlbiteRateFunction(reactionAlbite,system);
reactionAlbite.setRate(ReactionAlbiteRateFunction);

const Reaktoro::ReactionEquation equationQuartz("Quartz = SiO2@");
Reaktoro::Reaction reactionQuartz(equationQuartz, system);
reactionQuartz.setName("reactionQuartz");
Reaktoro::ReactionRateFunction ReactionQuartzRateFunction = theQuartzRateFunction(reactionQuartz,system);
reactionQuartz.setRate(ReactionQuartzRateFunction);

const Reaktoro::ReactionEquation equationKaolinite("Kaolinite + 6*H+ = 5*H2O@ + 2*Al+3 + 2*SiO2@");
Reaktoro::Reaction reactionKaolinite(equationKaolinite, system);
reactionKaolinite.setName("reactionKaolinite");
Reaktoro::ReactionRateFunction ReactionKaoliniteRateFunction = theKaoliniteRateFunction(reactionKaolinite,system);
reactionKaolinite.setRate(ReactionKaoliniteRateFunction);

const std::vector<Reaktoro::Reaction> reactions{reactionMicrocline,reactionQuartz,reactionKaolinite,reactionAlbite};