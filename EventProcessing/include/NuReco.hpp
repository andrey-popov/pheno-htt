#pragma once

#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TLorentzVector.h>

#include <utility>


/**
 * \class NuReco
 * \brief Reconstructs neutrino using constrants from masses of top quark and W boson and
 * missing pt
 * 
 * Reconstruction algorithm is based on [1]. The two mass constraints define an ellipse in the
 * space of neutrino's three-momentum. A unique solution on this ellipse is chosen by minimizing
 * Euclidian distance to measured missing pt.
 * [1] B.A. Betchart, R. Demina, A. Harel, Nucl.Instrum.Meth. A736 (2014) 169 [arXiv:1305.1878]
 * 
 * The algorithm will deliver no solution if the two mass constraints are not compatible with each
 * other.
 */
class NuReco
{
private:
    using Matrix = ROOT::Math::SMatrix<double, 3>;
    using Vector = ROOT::Math::SVector<double, 3>;
    
public:
    /**
     * \brief Constructor
     * 
     * Provided masses of top quark and W boson will be used in the constraints.
     */
    NuReco(double massTop, double massW);
    
public:
    /**
     * \brief Returns a measure of compatibility between obtained solution and input missing pt
     * 
     * The compatibility is evaluated as Euclidian distance between the two transverse momenta.
     */
    double GetCompatibility() const;
    
    /**
     * \brief Returns momentum of reconstructed neutrino
     * 
     * The momentum is zero if reconstruction has not succeeded.
     */
    TLorentzVector const &GetSolution() const;
    
    /// Reconstructs neutrino for given momenta of lepton and b quark and missing pt
    void Reconstruct(TLorentzVector const &p4Lep, TLorentzVector const &p4B,
      TLorentzVector const &p4Miss) const;
    
    /**
     * \brief Reports status of reconstruction
     * 
     * Status 0 indicates a successful reconstruction. Status 1 means that the two mass constraints
     * are not consistent with each other and no solution can be found.
     */
    unsigned RecoStatus() const;
    
    /// Changes masses of top quark and W boson used in the constraints
    void SetMasses(double massTop, double massW);
    
private:
    /**
     * \brief Finds the solution that is most compatible with measured missing pt
     * 
     * Minimization algorithm exploits also the first derivative of the loss function.
     */
    std::pair<double, double> Minimize() const;
    
    /**
     * \brief Constructs a rotation matrix
     * 
     * Implements rotation about an axis of the coordinate system by the given angle. Axes x, y, z
     * are denoted by numbers 0, 1, 2.
     */
    static Matrix Rotation(unsigned axis, double angle);
    
    /**
     * \brief Solves the two mass constraints
     * 
     * Technically, matrix H is constructed (Sec. 2.5 of the paper). Returns a boolean that
     * indicates whether the two constraints can be met simultaneously.
     */
    bool SolveMassConstraints(TLorentzVector const &p4Lep, TLorentzVector const &p4B) const;
    
private:
    /// (Squared) masses of top quark and W boson (GeV^2)
    double massTop2, massW2;
    
    /**
     * \brief Transformation from parameterized solution for neutrino's momentum to lab coordinates
     * 
     * The solutions are represented with vector T = (cos(t), sin(t), 1)', where t is the
     * parameter.
     */
    mutable Matrix H;
    
    /**
     * \brief Matrices to compute the distance between the solution and input missing pt and the
     * derivative of this distance with respect to the parameter
     * 
     * The squared distance and its derivative are given by T'XT and T'MT respectively.
     */
    mutable Matrix X, M;
    
    /// Indicates whether the reconstruction has been successful
    mutable bool recoSuccess;
    
    /// Distance between obtained solution for neutrino momentum and input missing pt
    mutable double distance;
    
    /// Reconstructed four-momentum of the neutrino
    mutable TLorentzVector p4Nu;
};
