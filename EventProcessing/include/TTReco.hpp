#pragma once

#include <AnalysisPlugin.hpp>

#include <DelphesReaderBase.hpp>
#include <NuReco.hpp>

#include <string>
#include <vector>


class LJetsSelection;

class TH1;
class TH2;


/**
 * \class TTReco
 * 
 * A plugin to perform jet assignment in a semileptonic ttbar event
 * 
 * This plugin identifies reconstructed jets corresponding to the four quarks in the final state
 * of tt -> blv bqq. Neutrino is reconstructed with a dedicated class for each considered choice
 * of the corresponding b-quark jet. All possible ways to assign jets to the four quarks, while
 * respecting the b-tagging decisions, are considered, and for each interpretation a likelihood is
 * evaluated, which combines the measure of consistency of reconstruction of neutrino and
 * reconstructed masses of hadronically decaying top quark and W boson. The interpretation with the
 * highest rank is accepted.
 */
class TTReco: public AnalysisPlugin
{
public:
    /// Jets to be identified in the final state of a ttbar system
    enum class DecayJet
    {
        bTopLep,   ///< Jet from semileptonically decaying top quark
        bTopHad,   ///< Jet from hadronization of b quark from hadronically decaying top quark
        q1TopHad,  ///< Leading light-flavour jet from hadronically decaying top quark
        q2TopHad   ///< Subleading light-flavour jet from hadronically decaying top quark
    };
    
public:
    /// Constructor from pointers to required plugins and a path to file that defines likelihood
    TTReco(DelphesReaderBase const *reader, LJetsSelection const *selector,
      std::string const &likelihoodFile);
    
public:
    /**
     * Returns jet corresponding to the given quark in the final state tt -> blv bqq
     * 
     * The behaviour is undefined if reconstruction has failed.
     */
    Jet const &GetJet(DecayJet type) const;
    
    /// Returns four-momentum of the charged lepton from the t->blv decay
    TLorentzVector const &GetLeptonP4() const;
    
    /// Returns four-momentum of reconstructed neutrino from the t->blv decay
    TLorentzVector const &GetNeutrinoP4() const;
    
    /**
     * Returns rank of the accepted interpretation of the current event
     * 
     * If reconstruction has failed, returns -infinity.
     */
    double GetRank() const;
    
    /**
     * Returns status code showing whether the reconstruction has succeeded
     * 
     * Code 0 indicates a successful reconstruction, any other points to a failure.
     */
    unsigned GetRecoStatus() const;
    
    /// Computes and returns four-momentum of reconstructed leptonically decaying top quark
    TLorentzVector GetTopLepP4() const;
    
    /// Computes and returns four-momentum of reconstructed hadronically decaying top quark
    TLorentzVector GetTopHadP4() const;
    
    /**
     *  Sets jet selection
     * 
     * Only jets satisfying this selection are tried as decay products of the top quarks.
     */
    void SetJetSelection(double minPt, double maxAbsEta = std::numeric_limits<double>::infinity());
    
private:
    /**
     * Performs reconstruction of the current event
     * 
     * Considers all possible ways to choose four reconstructed jets and match them to decay
     * products of a pair of top quarks. Only b-tagged jets are considered for matching to the
     * b quarks. Jets must also satisfy a selection on pt and |eta|. For each considered
     * interpretation the neutrino is reconstructed, and the rank of the interpretation is
     * evaluated. The rank is defined as the logarithm of the combined likelihood for the neutrino
     * distance and (mt, mW).
     * 
     * The event is rejected if the reconstruction fails.
     */
    virtual bool ProcessEvent() override;
    
private:
    /// Non-owning pointer to reader plugin
    DelphesReaderBase const *reader;
    
    /// Non-owning pointer to plugin that performs event selection
    LJetsSelection const *selector;
    
    /// Object that performs reconstruction of neutrino
    NuReco nuReco;
    
    /// Histograms that define likelihood function
    TH1 *likelihoodNuDist;
    TH2 *likelihoodMassHad;
    
    /// Kinematic selection for jets
    double minPt, maxAbsEta;
    
    /**
     * Indices of jets that pass the selection on pt and |eta|
     * 
     * This vector is only used in the method ProcessEvent but placed here to avoid reallocation of
     * memory for each event.
     */
    std::vector<unsigned> selectedJetIndices;
    
    /**
     * Status code indicating success of failure of reconstruction
     * 
     * Code 0 indicates a successful reconstruction. Other values denote failures and describe
     * their reasons.
     */
    unsigned recoStatus;
    
    /**
     * Rank of the best interpretation constructed so far
     * 
     * When reconstruction of a new event starts, this variable is reset to -infinity. After all
     * interpretations in an event have been considered, it contains the rank of the best
     * interpretation.
     */
    double highestRank;
    
    /**
     * Non-owning pointers to jets identified as decay products of the top quarks
     * 
     * Pointers refer to jets in the collection provided by the jet reader.
     */
    Jet const *bTopLep, *bTopHad, *q1TopHad, *q2TopHad;
};
