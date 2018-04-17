#include "PhysicsTools/NanoAODJMARTools/interface/Recluster.h"
#include "fastjet/PseudoJet.hh"

SoftDropWrapper::SoftDropWrapper(float beta, float zcut, float R, float ptmin) :
  fastjet::contrib::SoftDrop(beta,zcut), R_(R), ptmin_(ptmin)
{
}
 
SoftDropWrapper::SoftDropWrapper() :
  fastjet::contrib::SoftDrop(0.0,0.1), R_(0.8), ptmin_(0.0)
{
}

std::vector<fastjet::PseudoJet> SoftDropWrapper::result( std::vector<TLorentzVector> const & particles )
{
  fastjet::JetDefinition jet_def (fastjet::antikt_algorithm, R_);
  std::vector<fastjet::PseudoJet> pjets;
  for( auto const & particle : particles ) {
    pjets.emplace_back(particle.Px(), particle.Py(), particle.Pz(), particle.E());
  }
  fastjet::ClusterSequence cs( pjets, jet_def );
  std::vector<fastjet::PseudoJet> sorted = fastjet::sorted_by_pt( cs.inclusive_jets( ptmin_ ) );
  return this->operator()( sorted );
  
}
