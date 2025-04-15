/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "waveCurrent.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace waveTheories
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(waveCurrent, 0);
addToRunTimeSelectionTable(waveTheory, waveCurrent, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

waveCurrent::waveCurrent(
    const word &subDictName,
    const fvMesh &mesh_)
    : waveTheory(subDictName, mesh_),
      H_(readScalar(coeffDict_.lookup("height"))),
      h_(readScalar(coeffDict_.lookup("depth"))),
      omega_(readScalar(coeffDict_.lookup("omega"))),
      period_(2 * PI_ / omega_),
      phi_(readScalar(coeffDict_.lookup("phi"))),
      k_(vector(coeffDict_.lookup("waveNumber"))),
      k1_ = linearWaveNumber1(),
      U_(vector(coeffDict_.lookup("U"))),
      #omegac_(omega_ + (k_ & U_)),
      K_(mag(k_)),
      K1_(mag(k1_)),
      Tsoft_(coeffDict_.lookupOrDefault<scalar>("Tsoft", period_)),
      debug_(Switch(coeffDict_.lookup("debug")))
{
    checkWaveDirection(k_);

    if
    (
        H_/2.0 - 4.0*1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
        - 1.0/Foam::tanh(K_*h_)) < 0
    )
    {
        if (debug_)
        {
            WarningIn
            (
                "label waveCurrent::eta(point x, scalar time)"
            ) << endl << "The validity of stokes second order is violated."
            << endl << "a_1 < 4 a_2, being first and second order"
            << " amplitudes respectively." << endl << endl;
            Info << "a1 = " << H_/2.0 << " , a2 = "
                 << (1.0/16.0*K_*sqr(H_)*(3.0/Foam::pow(Foam::tanh(K_*h_),3.0)
                    - 1.0/Foam::tanh(K_*h_))) << endl;
        }
    }
}


void waveCurrent::printCoeffs()
{
    Info << "Loading wave theory: " << typeName << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


scalar waveCurrent::factor(const scalar& time) const
{
    scalar factor(1.0);
    if (Tsoft_ > 0.0)
    {
        factor = Foam::sin(2*PI_/(4.0*Tsoft_)*Foam::min(Tsoft_, time));
    }

    return factor;
}

scalar waveCurrent::eta
(
    const point& x,
    const scalar& time
) const
{
    scalar arg(omega_*time - (k1_ & x) + phi_);

    scalar eta = (
                     H_/2.0*Foam::cos(arg) // First order contribution.
                   // Second order contribution.
                   + 1.0/16.0*K1_*sqr(H_)
                   *(
                        3.0/Foam::pow(Foam::tanh(K1_*h_),3.0)
                      - 1.0/Foam::tanh(K1_*h_)
                    )
                   *Foam::cos(2.0*arg)
                 )*factor(time)  // Hot-starting.
                 + seaLevel_;      // Adding sea level.

    return eta;
}


//scalar waveCurrent::ddxPd
//(
//    const point& x,
//    const scalar& time,
//    const vector& unitVector
//) const
//{
//    scalar Z(returnZ(x));
//    scalar arg(omega_*time - (k_ & x) + phi_);
//
//    scalar ddxPd(0);
//
//    ddxPd = (
//                rhoWater_*mag(g_)*K_*H_/2.0*Foam::cosh(K_*(Z + h_))
//                /Foam::cosh(K_*h_)*Foam::sin(arg)
//                + 1/4*rhoWater_*mag(g_)*pow(K_,2)*pow(H_,2)
//                /Foam::sinh(2*K_*h_)
//                *( 3*Foam::cosh(2*K_*(Z + h_))/pow(Foam::sinh(K_*h_),2) - 1)
//                *Foam::sin(2*arg)
//            )*factor(time);
//
//    return ddxPd;
//}


scalar waveCurrent::pExcess
(
    const point& x,
    const scalar& time
) const
{
	scalar res = 0;

	// Get arguments and local coordinate system
    scalar Z(returnZ(x));
    scalar arg(omega_*time - (k1_ & x) + phi_);

    // First order contribution
    res = rhoWater_*mag(g_)*H_/2.0*Foam::cosh(K1_*(Z + h_))
        /Foam::cosh(K1_*h_)*Foam::cos(arg);

    // Second order contribution
    res += 1.0/8.0*rhoWater_*mag(g_)*K1_*Foam::sqr(H_)/Foam::sinh(2.0*K1_*h_)
        *((3.0*Foam::cosh(2.0*K1_*(Z + h_))/Foam::sqr(Foam::sinh(K1_*h_)) - 1.0)*Foam::cos(2.0*arg)
         - Foam::cosh(2*K1_*(Z + h_)) + 1);

    // Apply the ramping-factor
    res *= factor(time);
    res += referencePressure();

    return res;
}



vector waveCurrent::U
(
    const point& x,
    const scalar& time
) const
{
    scalar Z(returnZ(x));
    scalar cel(omega_/K1_);
    scalar arg(omega_*time - (k1_ & x) + phi_);

    // First order contribution
    scalar Uhorz = PI_*H_/period_ *
                   Foam::cosh(K1_*(Z + h_))/Foam::sinh(K1_*h_) *
                   Foam::cos(arg);

    // Second order contribution
    Uhorz += 3.0/16.0*cel*Foam::sqr(K1_*H_)*Foam::cosh(2*K1_*(Z + h_))
            /Foam::pow(Foam::sinh(K1_*h_),4.0)*Foam::cos(2*arg)
             - 1.0/8.0*mag(g_)*sqr(H_)/(cel*h_);

    // Current velocity in x dirction
    Uhorz += U_.x();

    // First order contribution
    scalar Uvert = - PI_*H_/period_ *
                   Foam::sinh(K1_*(Z + h_))/Foam::sinh(K1_*h_) *
                   Foam::sin(arg);

    // Second order contribution
    Uvert += - 3.0/16.0*cel*sqr(K1_*H_)*Foam::sinh(2*K1_*(Z + h_))
            /Foam::pow(Foam::sinh(K1_*h_), 4.0)*Foam::sin(2*arg);

    // Current velocity in y dirction
    Uvert += U_.y();

    // Multiply by the time stepping factor
    Uvert *= factor(time);
    Uhorz *= factor(time);

    // Generate the velocity vector
    // Note "-" because of "g" working in the opposite direction
    return Uhorz*k1_/K1_ - Uvert*direction_;
}

scalar waveCurrent::linearWaveNumber1()
{
    scalar lower(0.0);  // 设置区间的下界

    // 设置区间的上界
    scalar upper = Foam::max
        (
            4.0*PI_/( period_*Foam::sqrt( Foam::mag(G_)*depth_)),
            2.0*PI_/( Foam::pow( period_, 2.0))
        );

    scalar middle(0.5*(lower + upper));  // 计算中点

    scalar tanhMax(100);  // 最大值限制，用于计算tanh

    // 计算目标方程的初始值
    scalar valLower = Foam::pow(omega_ - lower * U_, 2.0) - lower * Foam::mag(G_) * Foam::tanh(lower * depth_);
    scalar valUpper = Foam::pow(omega_ - upper * U_, 2.0) - upper * Foam::mag(G_) * Foam::tanh(upper * depth_);
    scalar valMiddle = Foam::pow(omega_ - middle * U_, 2.0) - middle * Foam::mag(G_) * Foam::tanh(middle * depth_);

    // 使用二分法进行迭代
    while (true)
    {
        // 判断 valLower 和 valMiddle 的符号
        if (Foam::sign(valLower) == Foam::sign(valMiddle))
        {
            lower = middle;
            valLower = valMiddle;
        }
        else
        {
            upper = middle;
            valUpper = valMiddle;
        }

        // 计算新的 middle 值
        middle = 0.5 * (lower + upper);

        // 计算新的 valMiddle
        valMiddle = Foam::pow(omega_ - middle * U_, 2.0) - middle * Foam::mag(G_) * Foam::tanh(middle * depth_);

        // 检查收敛条件
        if
        (
            Foam::mag(valMiddle) < 1.0e-13 ||  // 当值足够小
            Foam::mag(valLower - valUpper) / middle < 1.0e-13  // 或者值变化足够小
        )
        {
            break;  // 满足收敛条件，跳出循环
        }
    }

    return middle;  // 返回计算得到的波数 k1
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace waveTheories
} // End namespace Foam

// ************************************************************************* //
