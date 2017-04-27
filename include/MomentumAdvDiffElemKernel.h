/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef MOMENTUMADVDIFFELEMKERNEL_H
#define MOMENTUMADVDIFFELEMKERNEL_H

#include "Kernel.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Kokkos_Core.hpp>

namespace sierra {
namespace nalu {

class SolutionOptions;
class MasterElement;
class ElemDataRequests;
class ScratchViews;

/** Advection diffusion term for momentum equation (velocity DOF)
 */
template<typename AlgTraits>
class MomentumAdvDiffElemKernel: public Kernel
{
public:
  MomentumAdvDiffElemKernel(
    const stk::mesh::BulkData&,
    SolutionOptions&,
    VectorFieldType*,
    ScalarFieldType*,
    ElemDataRequests&);

  virtual ~MomentumAdvDiffElemKernel();

  /** Execute the kernel within a Kokkos loop and populate the LHS and RHS for
   *  the linear solve
   */
  virtual void execute(
    SharedMemView<double**>&,
    SharedMemView<double*>&,
    stk::mesh::Entity,
    ScratchViews&);

private:
  MomentumAdvDiffElemKernel() = delete;

  VectorFieldType *velocityNp1_{nullptr};
  VectorFieldType *coordinates_{nullptr};
  ScalarFieldType *viscosity_{nullptr};
  GenericFieldType *massFlowRate_{nullptr};

  const int* lrscv_;

  const double includeDivU_;

  // fixed scratch space
  Kokkos::View<double[AlgTraits::numScsIp_][AlgTraits::nodesPerElement_]> v_shape_function_{"v_shape_function"};
  Kokkos::View<double[AlgTraits::nDim_]> v_uIp_{"v_uIp"};
};

}  // nalu
}  // sierra

#endif /* MOMENTUMADVDIFFELEMKERNEL_H */