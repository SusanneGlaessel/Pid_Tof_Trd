/** @file   ContainerTof.h
    @class  PidTof::Container
    @author Viktor Klochkov (klochkov44@gmail.com)
    @author Ilya Selyuzhenkov (ilya.selyuzhenkov@gmail.com)
    @brief  Class to store Tof PID information
*/

#ifndef CENTRALITY_CONTAINER_H
#define CENTRALITY_CONTAINER_H

#include <map>

typedef unsigned int uint;

namespace PidTof {

class ContainerTof {

 public:
  Container() = default;

  void SetBayesianProbability(std::map<uint, double>&& prob) { prob_ = prob; }
  void SetSigma(std::map<uint, double>&& sigma) { sigma_ = sigma; }

  double GetBayesianProbability(uint num) const {
    auto find = prob_.find(num);
    return find != prob_.end() ? find->second : -1;
  }

  double GetSigma(uint num) const {
    auto find = sigma_.find(num);
    return find != sigma_.end() ? find->second : -1;
  }

 private:
  std::map<uint, double> prob_;
  std::map<uint, double> sigma_;

  //     ClassDef(Container, 1);
};
}// namespace PidTof

#endif//CENTRALITY_CONTAINER_H
