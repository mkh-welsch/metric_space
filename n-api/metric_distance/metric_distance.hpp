/* Signal Empowering Technology   
                       
presents

███╗   ███╗███████╗████████╗██████╗ ██╗ ██████╗    ██████╗ ██╗███████╗████████╗ █████╗ ███╗   ██╗ ██████╗███████╗    
████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██║██╔════╝    ██╔══██╗██║██╔════╝╚══██╔══╝██╔══██╗████╗  ██║██╔════╝██╔════╝    
██╔████╔██║█████╗     ██║   ██████╔╝██║██║         ██║  ██║██║███████╗   ██║   ███████║██╔██╗ ██║██║     █████╗      
██║╚██╔╝██║██╔══╝     ██║   ██╔══██╗██║██║         ██║  ██║██║╚════██║   ██║   ██╔══██║██║╚██╗██║██║     ██╔══╝      
██║ ╚═╝ ██║███████╗   ██║   ██║  ██║██║╚██████╗    ██████╔╝██║███████║   ██║   ██║  ██║██║ ╚████║╚██████╗███████╗    
╚═╝     ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝ ╚═════╝    ╚═════╝ ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═══╝ ╚═════╝╚══════╝                                                                                                                                         
                                                                                          Michael Welsch (c) 2018.
                                                                                                   
a library for metrics / distance functions


This Source Code Form is subject to the terms of the PANDA GmbH
License. You are not allowed to use or edit the code without license.

*/

#ifndef _METRIC_DISTANCE_HPP
#define _METRIC_DISTANCE_HPP
#include <vector>
#include <limits>
#include "metric_factory.hpp"

// namespace metric
// {

// namespace distance
// {

// /*** Euclidian (L2) Metric ***/
// template <typename Container>
// struct Euclidian
// {
//   typename Container::value_type
//   operator()(const Container &a, const Container &b) const;
// };

// /***  Manhatten/Cityblock (L1) Metric ***/
// template <typename Container>
// struct Manhatten
// {
//   typename Container::value_type
//   operator()(const Container &a, const Container &b) const;
// };

// /*** Minkowski (L general) Metric ***/
// template <typename Container>
// struct P_norm
// {
//   typename Container::value_type
//   operator()(const Container &a, const Container &b) const;
//   typename Container::value_type p;
//   P_norm(typename Container::value_type p) : p(p){};
// };

// /*** Minkowski Metric (L... / P_Norm) ***/
// template <typename Container>
// struct Euclidian_thresholded
// {
//   typename Container::value_type
//   operator()(const Container &a, const Container &b) const;
//   typename Container::value_type thres;
//   typename Container::value_type factor;
//   Euclidian_thresholded(typename Container::value_type thres = 1000, typename Container::value_type factor = 3000) : thres(thres), factor(factor) {}
// };

// /*** Cosine Metric ***/
// template <typename Container>
// struct Cosine
// {
//   typename Container::value_type
//   operator()(const Container &a, const Container &b) const;
// };

// /*** structural similartiy (for images) ***/
// template <typename Container>
// struct SSIM
// {
//   double
//   operator()(const Container &img1, const Container &img2, double const &dynamic_range = 255, double const &masking = 2.0) const;
// };

// /*** Time Warp Elastic Distance (for curves) ***/
// template <typename Container>
// struct TWED
// {
//   typename Container::value_type
//   operator()(const Container &As, const Container &Bs, typename Container::value_type const &penalty = 0, typename Container::value_type const &elastic = 1) const;
// #ifdef _BLAZE_BLAZE_H_
//   typename Container::value_type
//   operator()(blaze::CompressedVector<Container::value_type> const &As, blaze::CompressedVector<Container::value_type> const &Bs, Container::value_type const &penalty = 0, Container::value_type const &elastic = 1, bool is_zero_padded = false) const;
// #endif
// };

// /*** structural similartiy (for images) ***/
// template <typename Container>
// struct EMD
// {
//   typename Container::value_type
//   operator()(const Container &Pc,
//              const Container &Qc,
//              const std::vector<std::vector<typename Container::value_type>> &C,
//              typename Container::value_type maxC = std::numeric_limits<typename Container::value_type>::min(),
//              typename Container::value_type extra_mass_penalty = -1,
//              std::vector<std::vector<typename Container::value_type>> *F = NULL) const;
// };

// /*** edit distance (for strings) ***/
// template <typename Container>
// struct Edit
// {
//     int
//     operator()(const Container &str1, const Container &str2) const;
// };

// } // namespace distance
// } // namespace metric
#include "metric_distance.cpp" // include the implementation

#endif //_METRIC_DISTANCE_HPP
