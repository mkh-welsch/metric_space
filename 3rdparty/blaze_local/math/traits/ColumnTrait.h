//=================================================================================================
/*!
//  \file blaze/math/traits/ColumnTrait.h
//  \brief Header file for the column trait
//
//  Copyright (C) 2012-2018 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_TRAITS_COLUMNTRAIT_H_
#define _BLAZE_MATH_TRAITS_COLUMNTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "../../util/InvalidType.h"
#include "../../util/mpl/If.h"
#include "../../util/mpl/Or.h"
#include "../../util/typetraits/Decay.h"
#include "../../util/typetraits/IsConst.h"
#include "../../util/typetraits/IsReference.h"
#include "../../util/typetraits/IsVolatile.h"


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the ColumnTrait class.
// \ingroup math_traits
//
// \section columntrait_general General
//
// The ColumnTrait class template offers the possibility to select the resulting data type when
// creating a view on a specific column of a dense or sparse matrix. ColumnTrait defines the nested
// type \a Type, which represents the resulting data type of the column operation. In case the
// given data type is not a dense or sparse matrix type, the resulting data type \a Type is
// set to \a INVALID_TYPE. Note that \a const and \a volatile qualifiers and reference modifiers
// are generally ignored.
//
//
// \section columntrait_specializations Creating custom specializations
//
// Per default, ColumnTrait supports all matrix types of the Blaze library (including views and
// adaptors). For all other data types it is possible to specialize the ColumnTrait template. The
// following example shows the according specialization for the DynamicMatrix class template:

   \code
   template< typename T1, bool SO, size_t... CCAs >
   struct ColumnTrait< DynamicMatrix<T1,SO>, CCAs... >
   {
      using Type = DynamicVector<T1,true>;
   };
   \endcode

// \n \section columntrait_examples Examples
//
// The following example demonstrates the use of the ColumnTrait template, where depending on
// the given matrix type the resulting column type is selected:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the column type of a column-major dynamic matrix
   using MatrixType1 = blaze::DynamicMatrix<int,columnMajor>;
   using ResultType1 = typename blaze::ColumnTrait<MatrixType1>::Type;

   // Definition of the column type for the 1st column of a row-major static matrix
   using MatrixType2 = blaze::StaticMatrix<int,3UL,4UL,rowMajor>;
   using ResultType2 = typename blaze::ColumnTrait<MatrixType2,1UL>::Type;
   \endcode
*/
template< typename MT       // Type of the matrix
        , size_t... CCAs >  // Compile time column arguments
struct ColumnTrait
{
 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   struct Failure { using Type = INVALID_TYPE; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename If_< Or< IsConst<MT>, IsVolatile<MT>, IsReference<MT> >
                            , ColumnTrait< Decay_<MT>, CCAs... >
                            , Failure >::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the ColumnTrait type trait.
// \ingroup math_traits
//
// The ColumnTrait_ alias declaration provides a convenient shortcut to access the nested
// \a Type of the ColumnTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename ColumnTrait<MT>::Type;
   using Type2 = ColumnTrait_<MT>;
   \endcode
*/
template< typename MT       // Type of the matrix
        , size_t... CCAs >  // Compile time column arguments
using ColumnTrait_ = typename ColumnTrait<MT,CCAs...>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
