/** @file   ConstantsTrd.h
    @author Susanne Glaessel (glaessel@ikf.uni-frankfurt.de)
    @brief  Some constants and enumerators
*/

#ifndef ConstantsTrd_H
#define ConstantsTrd_H

#include "TString.h"
#include <array>
#include <limits>

using std::make_pair;

namespace PidTrdParticles {
enum eParticle {
  kProton = 2212,
  kPionPos = 211,
  kKaonPos = 321,
  kDeutron = 1000010020,
  kTriton = 1000010030,
  kHe3 = 1000020030,
  kHe4 = 1000020040,
  kElectronPos = -11,
  kMuonPos = -13,
  kBgPos = 1,
  kAntiProton = -2212,
  kPionNeg = -211,
  kKaonNeg = -321,
  kDeutronNeg = -1000010020,
  kTritonNeg = -1000010030,
  kHe3Neg = -1000020030,
  kHe4Neg = -1000020040,
  kElectronNeg = 11,
  kMuonNeg = 13,
  kBgNeg = -1
};
};
constexpr int NbinsMax = 100;
constexpr int NumberOfTrdLayers = 4;
constexpr int NumberOfPidsTrd = 10;
constexpr int NumberOfTruncMode = 5;
constexpr int NumberOfProbMode = 2;

constexpr float BwMom = 0.05;
constexpr int NbinsMom = 60;
constexpr float BwdEdx = 0.1;
constexpr int NbinsdEdx = 35;

inline std::array<std::pair<Int_t, TString>, NumberOfPidsTrd> pid_codes_trd_{
    make_pair(2212, "p"),
    make_pair(211, "pi"),
    make_pair(321, "K"),
    make_pair(1000010020, "d"),
    make_pair(1000010030, "t"),
    make_pair(1000020030, "he3"),
    make_pair(1000020040, "he4"),
    make_pair(-11, "e"),
    make_pair(-13, "mu"),
    make_pair(1, "bg")};

inline TString dirname_tracks_ = "dEdx_tracks";
inline TString dirname_hits_ = "dEdx_hits";
inline std::array<TString, NumberOfTrdLayers> dirname_nhits_ = {
    "1hit",
    "2hits",
    "3hits",
    "4hits"};
inline std::array<TString, NumberOfTruncMode> histtitle_mode_ = {
    "truncation mode: 1 hit",
    "truncation mode: 2 hits",
    "truncation mode: 3 hits",
    "truncation mode: 4 hits"};

#endif// Constants_H
