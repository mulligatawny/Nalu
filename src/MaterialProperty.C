/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "MaterialProperty.h"

#include "Realm.h"
#include "MaterialPropertys.h"
#include "NaluEnv.h"

// props; algs, evaluators and data
#include "property_evaluator/MaterialPropertyData.h"
#include "property_evaluator/ReferencePropertyData.h"
#include "property_evaluator/PropertyEvaluator.h"

// yaml for parsing..
#include <yaml-cpp/yaml.h>

// basic c++
#include <stdexcept>
#include <map>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// MaterialProperty - material property
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
MaterialProperty::MaterialProperty(
  MaterialPropertys &materialPropertys, 
  const std::string materialBlockName)
  : materialPropertys_(materialPropertys),
    materialBlockName_(materialBlockName),
    propertyTableName_("na")
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
MaterialProperty::~MaterialProperty()
{   
  // delete prop data and evaluator
  std::map<PropertyIdentifier, MaterialPropertyData *>::iterator ii;
  for( ii=propertyDataMap_.begin(); ii!=propertyDataMap_.end(); ++ii )
    delete (*ii).second;
  
  std::map<PropertyIdentifier, PropertyEvaluator *>::iterator ie;
  for( ie=propertyEvalMap_.begin(); ie!=propertyEvalMap_.end(); ++ie )
    delete (*ie).second;
  
  // delete reference property map
  std::map<std::string, ReferencePropertyData *>::iterator irf;
  for( irf=referencePropertyDataMap_.begin(); irf!=referencePropertyDataMap_.end(); ++irf )
    delete (*irf).second;
}

//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
MaterialProperty::load(const YAML::Node &y_prop) 
{
  // extract the set of target names
  const YAML::Node targets = y_prop["target_name"];
  if (targets.Type() == YAML::NodeType::Scalar) {
    targetNames_.resize(1);
    targetNames_[0] = targets.as<std::string>() ;
  }
  else {
    targetNames_.resize(targets.size());
    for (size_t i=0; i < targets.size(); ++i) {
      targetNames_[i] = targets[i].as<std::string>() ;
    }
  }
  
  // has a table?
  if ( y_prop["table_file_name"] ) {
    propertyTableName_ = y_prop["table_file_name"].as<std::string>() ;
  }
  
  // property constants
  const YAML::Node y_constants = expect_map(y_prop, "constant_specification", true);
  if ( y_constants ) {
    universalConstantMap_ = y_constants.as<std::map<std::string, double> >() ;
  }
  
  // reference quantities
  const YAML::Node y_refs = expect_sequence(y_prop, "reference_quantities", true);
  if (y_refs) {
    for (size_t ispec = 0; ispec < y_refs.size(); ++ispec) {
      const YAML::Node y_ref = y_refs[ispec] ;
      
      // new the info object
      ReferencePropertyData *refData = new ReferencePropertyData();
      
      // extract the data; name and mw are always required
      get_required(y_ref, "species_name", refData->speciesName_);
      get_required(y_ref, "mw", refData->mw_);
      get_if_present(y_ref, "mass_fraction", refData->massFraction_);
      get_if_present(y_ref, "stoichiometry", refData->stoichiometry_);
      get_if_present(y_ref, "primary_mass_fraction", refData->primaryMassFraction_);
      get_if_present(y_ref, "secondary_mass_fraction", refData->secondaryMassFraction_);
      
      // store the map
      referencePropertyDataMap_[refData->speciesName_] = refData;      
    }
  }
  
  // extract the sequence of types
  const YAML::Node y_specs = expect_sequence(y_prop, "specifications", false);
  if (y_specs) {
    for (size_t ispec = 0; ispec < y_specs.size(); ++ispec) {
      const YAML::Node y_spec = y_specs[ispec];
      
      // new the info object
      MaterialPropertyData *matData = new MaterialPropertyData();
      
      std::string thePropName;
      std::string thePropType;
      get_required(y_spec, "name", thePropName);
      get_required(y_spec, "type", thePropType);
      
      // extract propery enum
      bool foundIt = false;
      PropertyIdentifier thePropEnum = PropertyIdentifier_END;
      for ( int k=0; k < PropertyIdentifier_END; ++k ) {
        if ( thePropName == PropertyIdentifierNames[k] ) {
          thePropEnum = PropertyIdentifier(k);
          foundIt = true;
          break;
        }
      }
      if ( !foundIt ) {
        throw std::runtime_error("could not find property name from enum list");
      }
      
      if ( thePropType == "constant" ) {
        matData->type_ = CONSTANT_MAT;
        
        // check for standard constant
        const YAML::Node standardConst = y_spec["value"];
        if ( standardConst ) {
          double theValue = 0.0;
          theValue = standardConst.as<double>() ;
          matData->constValue_ = theValue;
          NaluEnv::self().naluOutputP0() << thePropName
                                         << " is a constant property: " << theValue << std::endl;
        }
        else {
          // check for possible species Cp_k specification
          const YAML::Node coeffDeclare = y_spec["coefficient_declaration"];
          if (coeffDeclare) {
            const YAML::Node y_coeffds = expect_sequence(y_spec, "coefficient_declaration", true);
            if (y_coeffds) {
              for (size_t icoef = 0; icoef < y_coeffds.size(); ++icoef) {
                const YAML::Node y_coeffd = y_coeffds[icoef];
                std::string speciesName = "";
                get_required(y_coeffd, "species_name", speciesName);
                // standard coefficients single in size
                const YAML::Node coeffds = y_coeffd["coefficients"];
                if ( coeffds ) {
                  double theValue = 0.0;
                  theValue = coeffds.as<double>() ;
                  matData->cpConstMap_[speciesName] = theValue;
                }
                else {
                  throw std::runtime_error("no coefficients specified for Cp_k");
                }
                
                // standard hf single in size
                const YAML::Node hfs = y_coeffd["heat_of_formation"];
                if ( hfs ) {
                  double theValue = 0.0;
                  theValue = hfs.as<double>() ;
                  matData->hfConstMap_[speciesName] = theValue;
                }
                else {
                  NaluEnv::self().naluOutputP0() << "default heat of formation set to zero" << std::endl;
                  matData->hfConstMap_[speciesName] = 0.0;
                }
              }
            }
          }
          else {
            throw std::runtime_error("Cp must be either a constant single value or a set of constant coefficients provided");
          }
        }
      }
      else if ( thePropType == "mixture_fraction") {
        double primaryVal;
        double secondaryVal;
        get_required(y_spec, "primary_value", primaryVal);
        get_required(y_spec, "secondary_value", secondaryVal);
        matData->type_ = MIXFRAC_MAT;
        matData->primary_ = primaryVal;
        matData->secondary_ = secondaryVal;
        NaluEnv::self().naluOutputP0() << thePropName
                                       << " is a mix frac prop: "
                                       << primaryVal << " "
                                       << secondaryVal << std::endl;
      }
      else if ( thePropType == "volume_of_fluid") {
        double phaseOne;
        double phaseTwo;
        get_required(y_spec, "phase_one", phaseOne);
        get_required(y_spec, "phase_two", phaseTwo);
        matData->type_ = VOF_MAT;
        matData->phaseOne_ = phaseOne;
        matData->phaseTwo_ = phaseTwo;
        NaluEnv::self().naluOutputP0() << thePropName
                                       << " is a vof prop: "
                                       << phaseOne << " "
                                       << phaseTwo << std::endl;
      }
      else if ( thePropType == "polynomial" ) {
        matData->type_ = POLYNOMIAL_MAT;
        
        // check for coeff declaration
        const YAML::Node y_coeffds = expect_sequence(y_spec, "coefficient_declaration", true);
        if (y_coeffds) {
          for (size_t icoef = 0; icoef < y_coeffds.size(); ++icoef) {
            const YAML::Node y_coeffd = y_coeffds[icoef];
            std::string speciesName = "";
            get_required(y_coeffd, "species_name", speciesName);
            
            // standard coefficients
            const YAML::Node coeffds = y_coeffd["coefficients"];
            if ( coeffds ) {
              const size_t coeffSize = coeffds.size();
              std::vector<double> tmpPolyCoeffs(coeffSize);
              tmpPolyCoeffs = coeffds.as<std::vector<double> >();
              matData->polynomialCoeffsMap_[speciesName] = tmpPolyCoeffs;
            }
            
            // low temperature coefficient
            const YAML::Node lowcoeffds = y_coeffd["low_coefficients"];
            if ( lowcoeffds ) {
              const size_t coeffSize = lowcoeffds.size();
              std::vector<double> tmpPolyCoeffs(coeffSize);
              tmpPolyCoeffs = lowcoeffds.as<std::vector<double> >();
              matData->lowPolynomialCoeffsMap_[speciesName] = tmpPolyCoeffs;
            }
	    
            // high temperature coefficient
            const YAML::Node highcoeffds = y_coeffd["high_coefficients"];
            if ( highcoeffds ) {
              const size_t coeffSize = highcoeffds.size();
              std::vector<double> tmpPolyCoeffs(coeffSize);
              tmpPolyCoeffs = highcoeffds.as<std::vector<double> >();
              matData->highPolynomialCoeffsMap_[speciesName] = tmpPolyCoeffs;
            }
          }
        }
        else  {
          // complain if both are null
          throw std::runtime_error("polynominal property did not provide a set of coefficients");
        }
      }
      else if ( thePropType == "ideal_gas_t" || thePropType == "ideal_gas_t_p" || thePropType == "ideal_gas_yk"  ) {
        // warning about deprecation (backward compatibility)
        NaluEnv::self().naluOutputP0() << "MaterialProperty::warning(): usage of: " << thePropType << " for: " << thePropName 
                                       << " has been deprecated " << std::endl
                                       << " The new ideal gas name is simply: ideal_gas" << std::endl;
        NaluEnv::self().naluOutputP0() << thePropName
                                       << " will be set as an ideal gas property (function of tranported independent variables) " << std::endl;
        matData->type_ = IDEAL_GAS_MAT;
      }
      else if ( thePropType == "ideal_gas" ) {
        NaluEnv::self().naluOutputP0() << thePropName
                                       << " is an ideal gas property (function of tranported independent variables) " << std::endl;      
        matData->type_ = IDEAL_GAS_MAT;
      }
      else if ( thePropType == "geometric" ) {
        matData->type_ = GEOMETRIC_MAT;
        NaluEnv::self().naluOutputP0() << thePropName
                                       << " is a geometric property: " << std::endl;
      }
      else if ( thePropType == "hdf5table") {
        matData->type_ = HDF5_TABLE_MAT;
        NaluEnv::self().naluOutputP0() << thePropName << " is a hdf5 table look-up property: " << std::endl;
 	
        // what we might expect 
        std::string tablePropName = "na"; 
        std::string auxVarName = "na"; std::string tableAuxVarName = "na";
        
        // extract possible vector
        const YAML::Node names = y_spec["independent_variable_set"];
        if (names.Type() == YAML::NodeType::Scalar) {
          matData->indVarName_.resize(1);
          matData->indVarName_[0] = names.as<std::string>() ;
        }
        else {
          matData->indVarName_.resize(names.size());
          for (size_t i=0; i < names.size(); ++i) {
            matData->indVarName_[i] = names[i].as<std::string>();
          }
        }
	
        // extract possible vector of independent table names
        const YAML::Node &tableNames = y_spec["table_name_for_independent_variable_set"];
        if (tableNames.Type() == YAML::NodeType::Scalar) {
          matData->indVarTableName_.resize(1);
          matData->indVarTableName_[0] = tableNames.as<std::string>();
        }
        else {
          matData->indVarTableName_.resize(names.size());
          for (size_t i=0; i < tableNames.size(); ++i) {
            matData->indVarTableName_[i] = tableNames[i].as<std::string>() ;
          }
        }
	
        // the following are all single in size
        get_if_present_no_default(y_spec, "table_name_for_property", tablePropName);
        get_if_present_no_default(y_spec, "aux_variables", auxVarName);
        get_if_present_no_default(y_spec, "table_name_for_aux_variables", tableAuxVarName);
	
        // set matData
        matData->auxVarName_ = auxVarName;
        matData->tablePropName_ = tablePropName;
        matData->tableAuxVarName_ = tableAuxVarName;
      }
      else if ( thePropType == "generic" ) {
        matData->type_ = GENERIC;
        get_required(y_spec, "generic_property_evaluator_name", matData->genericPropertyEvaluatorName_);
      }
      else {
        throw std::runtime_error("unknown property type");  
      }
      
      // store the map
      propertyDataMap_[thePropEnum] = matData;
    }  
  } 
}

//--------------------------------------------------------------------------
//-------- extract_universal_constant --------------------------------------
//--------------------------------------------------------------------------
void
MaterialProperty::extract_universal_constant( 
  const std::string name, double &value, const bool useDefault)
{
  std::map<std::string, double >::iterator it
    = universalConstantMap_.find(name);
  if ( it == universalConstantMap_.end() ) {
    // not found
    if ( useDefault ) {
      NaluEnv::self().naluOutputP0() << "WARNING: Reference value for " << name << " not found " << " using " << value << std::endl;
    }
    else {
      throw std::runtime_error("Realm::setup_property: reference value not found: " + name);
    }
  }
  else {
    value = (*it).second;
  }
}

} // namespace nalu
} // namespace Sierra
