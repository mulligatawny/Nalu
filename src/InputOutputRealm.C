/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "InputOutputRealm.h"
#include "Realm.h"
#include "SolutionOptions.h"

// transfer
#include "xfer/Transfer.h"

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Field.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <Ioss_SubSystem.h>

// c++
#include <string>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// InputOutputRealm - holder of data to be transferred to/from
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
InputOutputRealm::InputOutputRealm(Realms& realms, const YAML::Node & node)
  : Realm(realms, node)
{
  // nothing now
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
InputOutputRealm::~InputOutputRealm()
{
  for ( size_t k = 0; k < inputOutputFieldInfo_.size(); ++k )
    delete inputOutputFieldInfo_[k];
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
InputOutputRealm::initialize()
{
  // bare minimum to register fields and to extract from possible mesh file
  setup_post_processing_algorithms();
  register_io_fields();
  ioBroker_->populate_mesh();
  ioBroker_->populate_field_data();
  create_output_mesh();
  input_variables_from_mesh();
  initialize_post_processing_algorithms();
}

//--------------------------------------------------------------------------
//-------- register_io_fields ------------------------------------------------------
//--------------------------------------------------------------------------
void
InputOutputRealm::register_io_fields() {
  // register fields; extract vector of field/part; only nodal for now
  std::string velocityName = "velocity";
  for ( size_t k = 0; k < inputOutputFieldInfo_.size(); ++ k ) {
    const std::string fieldName = inputOutputFieldInfo_[k]->fieldName_;
    const int fieldSize = inputOutputFieldInfo_[k]->fieldSize_;
    const std::string fieldType = inputOutputFieldInfo_[k]->fieldType_;
    std::vector<std::string> targetNames = inputOutputFieldInfo_[k]->targetNames_;

    // sanity check on the type
    if ( fieldType != "node_rank" ) {
      throw std::runtime_error("Input/output realm only supports nodal_rank types");
    }

    // loop over target parts and declare/put the field
    for ( size_t j = 0; j < targetNames.size(); ++j ) {
      const std::string targetName = targetNames[j];
      stk::mesh::Part *targetPart = metaData_->get_part(targetName);
      if ( NULL == targetPart ) {
        throw std::runtime_error("Sorry, no part name found by the name: " + targetName + " for field: " + fieldName);
      }
      else {
        if ( fieldName.find(velocityName) != std::string::npos ) { //FIXME: require FieldType?
          VectorFieldType *velocity = &(metaData_->declare_field<VectorFieldType>(stk::topology::NODE_RANK, fieldName));
          stk::mesh::put_field_on_mesh(*velocity, *targetPart, fieldSize, nullptr);
        }
        else {
          stk::mesh::Field<double, stk::mesh::SimpleArrayTag> *theField
            = &(metaData_->declare_field< stk::mesh::Field<double, stk::mesh::SimpleArrayTag> >(stk::topology::NODE_RANK, fieldName));
          stk::mesh::put_field_on_mesh(*theField,*targetPart,fieldSize, nullptr);
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
InputOutputRealm::load(const YAML::Node & node)
{
  // call base class
  Realm::load(node);

  // now proceed with specific line commands to IO Realm
  const YAML::Node y_field = node["field_registration"];
  if (y_field) {

    // extract the sequence of types
    const YAML::Node y_specs = expect_sequence(y_field, "specifications", false);
    if (y_specs) {
      for (size_t ispec = 0; ispec < y_specs.size(); ++ispec) {
        const YAML::Node y_spec = y_specs[ispec] ;

        // find the name, size and type
        const YAML::Node fieldNameNode = y_spec["field_name"];
        const YAML::Node fieldSizeNode = y_spec["field_size"];
        const YAML::Node fieldTypeNode = y_spec["field_type"];

        if ( ! fieldNameNode )
          throw std::runtime_error("Sorry, field name must be provided");

        if ( ! fieldSizeNode )
          throw std::runtime_error("Sorry, field size must be provided");

        if ( ! fieldTypeNode )
          throw std::runtime_error("Sorry, field type must be provided");

        // new the data
        InputOutputInfo *theInfo = new InputOutputInfo();

        // push data to the info object
        theInfo->fieldName_ = fieldNameNode.as<std::string>() ;
        theInfo->fieldSize_ = fieldSizeNode.as<int>() ;
        theInfo->fieldType_ = fieldTypeNode.as<std::string>() ;

        const YAML::Node &targets = y_spec["target_name"];
        if (targets.Type() == YAML::NodeType::Scalar) {
          theInfo->targetNames_.resize(1);
          theInfo->targetNames_[0] = targets.as<std::string>() ;
        }
        else {
          theInfo->targetNames_.resize(targets.size());
          for (size_t i=0; i < targets.size(); ++i) {
            theInfo->targetNames_[i] = targets[i].as<std::string>() ;
          }
        }

        inputOutputFieldInfo_.push_back(theInfo);
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- populate_restart ------------------------------------------------
//--------------------------------------------------------------------------
double
InputOutputRealm::populate_restart(
  double &timeStepNm1, int &timeStepCount)
{
  return get_current_time();
}

//--------------------------------------------------------------------------
//-------- populate_external_variables_from_input --------------------------
//--------------------------------------------------------------------------
void
InputOutputRealm::populate_external_variables_from_input(
  const double currentTime)
{
  // only works for external field realm
  if ( type_ == "external_field_provider" && solutionOptions_->inputVarFromFileMap_.size() > 0 ) {
    std::vector<stk::io::MeshField> missingFields;
    const double foundTime = ioBroker_->read_defined_input_fields(currentTime, &missingFields);
    if ( missingFields.size() > 0 ) {
      for ( size_t k = 0; k < missingFields.size(); ++k) {
        NaluEnv::self().naluOutputP0() << "WARNING: Realm::populate_external_variables_from_input for field "
            << missingFields[k].field()->name()
            << " is missing; will default to IC specification" << std::endl;
      }
    }
    NaluEnv::self().naluOutputP0() << "Realm::populate_external_variables_from_input() candidate input time: "
                                   << foundTime << " for Realm: " << name() << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- compute_minimum_distance_to_wall  -------------------------------
//--------------------------------------------------------------------------
void InputOutputRealm::compute_wall_distance(const YAML::Node& wdist) {

    stk::mesh::PartVector fluid_parts_;
    stk::mesh::PartVector wall_parts_;
    std::string wall_dist_name_ = "NULL";

    auto fluid_partnames = wdist["fluid_parts"].as<std::vector<std::string>>();
    auto wall_partnames = wdist["wall_parts"].as<std::vector<std::string>>();

    if(wdist["wall_dist_name"]) {
        wall_dist_name_ = wdist["wall_dist_name"].as<std::string>();
    }

    fluid_parts_.resize(fluid_partnames.size());
    wall_parts_.resize(wall_partnames.size());

    for(size_t i=0; i<fluid_partnames.size(); i++) {
        stk::mesh::Part* part = meta_data.get_part(fluid_partnames[i]);
        fluid_parts_[i] = part;
    }

    for(size_t i=0; i<wall_partnames.size(); i++) {
        stk::mesh::Part* part = meta_data.get_part(wall_partnames[i]);
        wall_parts_[i] = part;
    }

    // Mesh meta and bulk data
    stk::mesh::BulkData & bulk_data = realm_.bulk_data();
    stk::mesh::MetaData & meta_data = realm_.meta_data();

    const unsigned nDim = meta_data.spatial_dimension();

    // Register fields and put on part
    VectorFieldType* coords = meta_data.get_field<VectorFieldType>(
            stk::topology::NODE_RANK, "coordinates");

    ScalarFieldType& ndtw = meta_data.declare_field<ScalarFieldType>(
            stk::topology::NODE_RANK, wall_dist_name_);


    for(auto part: fluid_parts_) {
        stk::mesh::put_field_on_mesh(*coords, *part, nDim, nullptr);
        stk::mesh::put_field_on_mesh(ndtw, *part, nullptr);
    }

    stk::mesh::Selector fluid_union = stk::mesh::selectUnion(fluid_parts_);
    stk::mesh::Selector wall_union = stk::mesh::selectUnion(wall_parts_);

    const stk::mesh::BucketVector& fluid_bkts = bulk_data.get_buckets(
            stk::topology::NODE_RANK, fluid_union);
    const stk::mesh::BucketVector& wall_bkts = bulk_data.get_buckets(
            stk::topology::NODE_RANK, wall_union);

    VectorFieldType* coords = meta_data.get_field<VectorFieldType>(
            stk::topology::NODE_RANK, "coordinates");
    ScalarFieldType* ndtw = meta_data.get_field<ScalarFieldType>(
            stk::topology::NODE_RANK, wall_dist_name_);


    // loop to compute ndtw
    for(size_t ib=0; ib < fluid_bkts.size(); ib++) {
        stk::mesh::Bucket& fbkt = *fluid_bkts[ib];
        double* xyz = stk::mesh::field_data(*coords, fbkt);
        double* wdist = stk::mesh::field_data(*ndtw, fbkt);

        for(size_t in=0; in < fbkt.size(); in++) {
            double min_dist = 1.0e16;
            for(size_t jb=0; jb < wall_bkts.size(); jb++) {
                stk::mesh::Bucket& wbkt = *wall_bkts[jb];
                double* wxyz = stk::mesh::field_data(*coords, wbkt);

                for(size_t jn=0; jn< wbkt.size(); jn++) {
                    double dist_calc = 0.0;
                    for(int j=0; j<nDim; j++) {
                        double dst = xyz[in*nDim+j] - wxyz[jn*nDim+j];
                        dist_calc += dst * dst;
                    }
                    if (dist_calc < min_dist) min_dist = dist_calc;
                }
            }
            wdist[in] = std::sqrt(min_dist);
        }
    }
}



} // namespace nalu
} // namespace Sierra
