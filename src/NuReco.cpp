#include <NuReco.hpp>

#include <Math/GenVector/VectorUtil.h>

#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>


NuReco::NuReco(double massTop, double massW):
    recoSuccess(false)
{
    SetMasses(massTop, massW);
}


double NuReco::GetCompatibility() const
{
    return distance;
}


TLorentzVector const &NuReco::GetSolution() const
{
    return p4Nu;
}


void NuReco::Reconstruct(TLorentzVector const &p4Lep, TLorentzVector const &p4B,
  TLorentzVector const &p4Miss) const
{
    // Construct matrix H, which provides solutions that meet the mass constraints. Abort
    //reconstruction in case of failure.
    bool const solvable = SolveMassConstraints(p4Lep, p4B);
    
    if (not solvable)
    {
        recoSuccess = false;
        p4Nu.SetXYZM(0, 0, 0, 0);
        return;
    }
    
    
    // Neutrino momentum in the lab coordinates can be found by applying transformation H to a
    //vector T = (cos(t), sin(t), 1)', where t parameterizes solutions. Now only need to find the
    //solution that matches measured missing pt best.
    
    
    // Missing pt in matrix form (V0 from Sec. 2.6.1 of the paper, with unconstrained z component)
    Matrix V0;  // initialized with zeros
    V0(0, 2) = p4Miss.Px();
    V0(1, 2) = p4Miss.Py();
    
    // Inverted error matrix for missing pt. Use simple identity matrix as the error matrix.
    Matrix Sigma2 = ROOT::Math::SMatrixIdentity();
    Sigma2(2, 2) = 0.;
    
    // Construct matrix X from Sec. 2.6.1. Compatibility between a solution for the neutrino's
    //momentum HT and the measured missing pt is quantified by the form T'XT. With the trivial
    //error matrix used here, this quantity coincides with squared Euclidian distance between the
    //two transverse momenta.
    Matrix const Lambda = V0 - H;
    X = ROOT::Math::Transpose(Lambda) * Sigma2 * Lambda;
    
    
    // Construct matrix M from Sec. 2.6.1. It is used in minimization of T'XT as the derivative of
    //this quantity with respect to parameter t is T'MT.
    Matrix D;  // initialized with zeros
    D(0, 1) = -1.;
    D(1, 0) = 1.;
    
    M = X * D;
    M += ROOT::Math::Transpose(M);
    
    
    // Find value of parameter t that gives the best solution
    auto const minimum = Minimize();
    
    
    // Translate the solution into lab coordinates and store additional information about it
    recoSuccess = true;
    
    Vector T;
    T[0] = std::cos(minimum.first);
    T[1] = std::sin(minimum.first);
    T[2] = 1.;
    
    Vector const p3Nu = H * T;
    p4Nu.SetXYZM(p3Nu[0], p3Nu[1], p3Nu[2], 0.);
    
    if (minimum.second >= 0.)
        distance = std::sqrt(minimum.second);
    else
    {
        // On rare occasions when squared distance is very close to zero, it might be negative due
        //to rounding errors. Preserve the sign to notify the user.
        distance = -std::sqrt(-minimum.second);
    }
}


unsigned NuReco::RecoStatus() const
{
    return (recoSuccess) ? 0 : 1;
}


void NuReco::SetMasses(double massTop, double massW)
{
    massTop2 = std::pow(massTop, 2);
    massW2 = std::pow(massW, 2);
}


std::pair<double, double> NuReco::Minimize() const
{
    Vector T;
    T[2] = 1.;  // this element will not change
    
    
    // The loss function can have up to two local minima. First will perform a rough scan over an
    //equidistant grid of points and identify consequitive pairs in which the derivative of the
    //loss function changes sign from negative to positive. They will provide approximate positions
    //of minima.
    
    unsigned const nPoints = 100;
    double const step = 2 * M_PI / nPoints;
    
    // Compute the derivative at the initial point
    double t = 0.;
    T[0] = std::cos(t);
    T[1] = std::sin(t);
    
    double prevDerivative = ROOT::Math::Similarity(M, T);
    
    // Positions of centres of identified pairs of points will be stored in an array. Do not use a
    //vector in order to avoid dynamic memory allocation. This is possible since the loss function
    //is known to have not more than two minima.
    std::array<double, 2> approxMinima;
    unsigned nMinima = 0;
    
    // Do the scan over grid. It wraps including the point at 2 pi to consider all consequitive
    //pairs.
    t = step;
    
    for (unsigned iPoint = 1; iPoint <= nPoints; ++iPoint, t += step)
    {
        // Compute derivative at the new point
        T[0] = std::cos(t);
        T[1] = std::sin(t);
        
        double const newDerivative = ROOT::Math::Similarity(M, T);
        
        
        // Check if signs of the derivatives are compatible with the hypothesis that there is a
        //minimum between the current and previous points
        if (prevDerivative < 0. and newDerivative > 0.)
        {
            // Save the position of the center of the pair
            approxMinima[nMinima] = t - step / 2;
            ++nMinima;
            
            if (nMinima == 2)
                break;
        }
        
        prevDerivative = newDerivative;
    }
    
    if (nMinima == 0)
        throw std::logic_error("NuReco::Minimize: Failed to find approximate position of a "
          "minimum in a grid search.");
    
    
    // Check all the approxate local minima and find position of the global one
    double minF = std::numeric_limits<double>::infinity();
    double tBest = 0.;
    
    for (unsigned i = 0; i < nMinima; ++i)
    {
        // Refine position of the current local minimum with a bisection search for the zero of the
        //derivative T'MT
        double t0 = approxMinima[i];
        double tMin = t0 - step / 2, tMax = t0 + step / 2;
        
        while (tMax - tMin > 1e-8)
        {
            t0 = (tMin + tMax) / 2;
            
            T[0] = std::cos(t0);
            T[1] = std::sin(t0);
            
            double const fD = ROOT::Math::Similarity(M, T);
            
            if (fD > 0.)
            {
                // If function is growing at t0, the minimum must be to the left of it
                tMax = t0;
            }
            else
                tMin = t0;
        }
        
        
        // Compare the value of the loss function at this local minimum with the current best
        //estimate for the global minimum
        t0 = (tMin + tMax) / 2;
        
        T[0] = std::cos(t0);
        T[1] = std::sin(t0);
        
        double const f = ROOT::Math::Similarity(X, T);
        
        if (f < minF)
        {
            minF = f;
            tBest = t0;
        }
    }
    
    return std::make_pair(tBest, minF);
}


NuReco::Matrix NuReco::Rotation(unsigned axis, double angle)
{
    Matrix r = ROOT::Math::SMatrixIdentity();
    
    double const c = std::cos(angle);
    double const s = std::sin(angle);
    
    switch (axis)
    {
        case 0:  // x
            r(1, 1) = c;
            r(1, 2) = -s;
            r(2, 1) = s;
            r(2, 2) = c;
            break;
        
        case 1:  // y
            r(0, 0) = c;
            r(0, 2) = s;
            r(2, 0) = -s;
            r(2, 2) = c;
            break;
        
        case 2:  // z
            r(0, 0) = c;
            r(0, 1) = -s;
            r(1, 0) = s;
            r(1, 1) = c;
            break;
        
        default:
        {
            std::ostringstream message;
            message << "NuReco::Rotation: Illegal index of rotation axis (" << axis << ").";
            throw std::runtime_error(message.str());
        }
    }
    
    return r;
}


bool NuReco::SolveMassConstraints(TLorentzVector const &p4Lep, TLorentzVector const &p4B) const
{
    // Implementation closely follows the paper
    
    // Cosine and sine of the angle between three-momenta of lepton and b quark (Sec. 2.3)
    double const c = ROOT::Math::VectorUtil::CosTheta(p4Lep, p4B);
    double const s = std::sqrt(1 - c * c);
    
    // Eq. 2 from the paper
    double const x0p = -(massTop2 - massW2 - p4B.M2()) / (2 * p4B.E());
    double const x0 = -(massW2 - p4Lep.M2()) / (2 * p4Lep.E());
    
    // Relativistic beta
    double const betaLep = p4Lep.Beta(), betaB = p4B.Beta();
    
    // $\epsilon^2$ from the end of Sec. 2.2
    double const epsilon2 = massW2 * (1 - betaLep * betaLep);
    
    // Eq. 6 and expression for Sy from Sec. 2.4.2
    double const Sx = (x0 * betaLep - p4Lep.P() * (1 - betaLep * betaLep)) / (betaLep * betaLep);
    double const Sy = (x0p / betaB - c * Sx) / s;
    
    // $\omega$ from Sec. 2.4.2 (the variable with minus sign is not used)
    double const omega = (betaLep / betaB - c) / s;
    
    // Expressions from Sec. 2.4.2
    double const Omega2 = omega * omega + (1 - betaLep * betaLep);
    double const x1 = Sx - (Sx + omega * Sy) / Omega2;
    double const y1 = Sy - (Sx + omega * Sy) * omega / Omega2;
    double const Z2 = x1 * x1 * Omega2 - std::pow(Sy - omega * Sx, 2) -
      (massW2 - x0 * x0 - epsilon2);
    
    
    if (Z2 < 0.)
    {
        // The two mass constraints are not consistent with each other
        return false;
    }
    
    double const Z = std::sqrt(Z2);
    
    
    // $\tilde H$ matrix from Sec. 2.4.2
    Matrix Htilde;  // initialized with zeros
    
    Htilde(0, 0) = Z / std::sqrt(Omega2);
    Htilde(0, 2) = x1 - p4Lep.P();
    Htilde(1, 0) = omega * Z / std::sqrt(Omega2);
    Htilde(1, 2) = y1;
    Htilde(2, 1) = Z;
    
    
    // Construct rotation matrix R from Sec. 2.5. Implementation follows the reference Python
    //implementation from the Appendix.
    Matrix const Rz = Rotation(2, -p4Lep.Phi());
    Matrix const Ry = Rotation(1, M_PI / 2 - p4Lep.Theta());
    
    Vector p3B;
    p3B[0] = p4B.Px();
    p3B[1] = p4B.Py();
    p3B[2] = p4B.Pz();
    
    Vector const p3BRotated = Ry * Rz * p3B;
    Matrix const Rx = Rotation(0, -std::atan2(p3BRotated[2], p3BRotated[1]));
    
    Matrix const R = ROOT::Math::Transpose(Rx * Ry * Rz);
    
    
    // Finally, build matrix H from Sec. 2.5
    H = R * Htilde;
    
    return true;
}
