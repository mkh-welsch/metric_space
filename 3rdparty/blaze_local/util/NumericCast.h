//=================================================================================================
/*!
//  \file blaze/util/NumericCast.h
//  \brief Cast operators for numeric types
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

#ifndef _BLAZE_UTIL_NUMERICCAST_H_
#define _BLAZE_UTIL_NUMERICCAST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <limits>
#include "../util/Assert.h"
#include "../util/constraints/Numeric.h"
#include "../util/EnableIf.h"
#include "../util/Exception.h"
#include "../util/FalseType.h"
#include "../util/mpl/And.h"
#include "../util/mpl/Bool.h"
#include "../util/mpl/Not.h"
#include "../util/mpl/Or.h"
#include "../util/TrueType.h"
#include "../util/typetraits/IsFloatingPoint.h"
#include "../util/typetraits/IsIntegral.h"
#include "../util/typetraits/IsSame.h"
#include "../util/typetraits/IsSigned.h"
#include "../util/typetraits/IsUnsigned.h"


namespace blaze {

//=================================================================================================
//
//  AUXILIARY CLASS DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary type trait for the numeric_cast() function template.
// \ingroup util
*/
template< typename To, typename From >
struct IsCriticalIntIntConversion
{
 public:
   //**********************************************************************************************
   enum : size_t { value = And< IsIntegral<To>
                              , IsIntegral<From>
                              , Not< IsSame<To,From> > >::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary type trait for the numeric_cast() function template.
// \ingroup util
*/
template< typename To, typename From >
struct IsCriticalFloatIntConversion
{
 public:
   //**********************************************************************************************
   enum : size_t { value = And< IsIntegral<To>
                              , IsFloatingPoint<From> >::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary type trait for the numeric_cast() function template.
// \ingroup util
*/
template< typename To, typename From >
struct IsCriticalFloatFloatConversion
{
 public:
   //**********************************************************************************************
   enum : size_t { value = And< IsFloatingPoint<To>
                              , IsFloatingPoint<From>
                              , Bool< ( sizeof(To) < sizeof(From) ) > >::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary type trait for the numeric_cast() function template.
// \ingroup util
*/
template< typename To, typename From >
struct IsUncriticalConversion
{
 public:
   //**********************************************************************************************
   enum : size_t { value = !IsCriticalIntIntConversion<To,From>::value   &&
                           !IsCriticalFloatIntConversion<To,From>::value &&
                           !IsCriticalFloatFloatConversion<To,From>::value };
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC CAST OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Numeric cast operators */
//@{
template< typename To, typename From > inline To numeric_cast( From from );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend for numeric_cast() for uncritical conversions.
// \ingroup util
//
// \param from The numeric value to be converted.
// \return The converted value.
*/
template< typename To, typename From >
inline EnableIf_< IsUncriticalConversion<To,From>, To >
   numeric_cast_backend( From from ) noexcept
{
   return from;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of numeric_cast() for critical conversions between integral types.
// \ingroup util
//
// \param from The numeric value to be converted.
// \return The converted value.
// \exception std::overflow_error Invalid numeric cast (overflow).
// \exception std::underflow_error Invalid numeric cast (underflow).
*/
template< typename To, typename From >
inline EnableIf_< IsCriticalIntIntConversion<To,From>, To >
   numeric_cast_backend( From from )
{
   if( ( sizeof(To) < sizeof(From) || ( IsSigned<To>::value && IsUnsigned<From>::value ) ) &&
       ( from > From( std::numeric_limits<To>::max() ) ) )
      BLAZE_THROW_OVERFLOW_ERROR( "Invalid numeric cast (overflow)" );

   if( IsSigned<From>::value &&
       ( from < From( std::numeric_limits<To>::min() ) ) )
      BLAZE_THROW_UNDERFLOW_ERROR( "Invalid numeric cast (underflow)" );

   BLAZE_INTERNAL_ASSERT( from == From( To( from ) ), "Numeric cast failed" );

   return from;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of numeric_cast() for critical conversions from floating point to integral types.
// \ingroup util
//
// \param from The numeric value to be converted.
// \return The converted value.
// \exception std::overflow_error Invalid numeric cast (overflow).
// \exception std::underflow_error Invalid numeric cast (underflow).
*/
template< typename To, typename From >
inline EnableIf_< IsCriticalFloatIntConversion<To,From>, To >
   numeric_cast_backend( From from )
{
   using std::trunc;

   if( from > From( std::numeric_limits<To>::max() ) )
      BLAZE_THROW_OVERFLOW_ERROR( "Invalid numeric cast (overflow)" );

   if( from < From( std::numeric_limits<To>::min() ) )
      BLAZE_THROW_UNDERFLOW_ERROR( "Invalid numeric cast (underflow)" );

   BLAZE_INTERNAL_ASSERT( trunc( from ) == From( To( from ) ), "Numeric cast failed" );

   return from;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of numeric_cast() for critical conversions between floating point types.
// \ingroup util
//
// \param from The numeric value to be converted.
// \return The converted value.
// \exception std::overflow_error Invalid numeric cast (overflow).
// \exception std::underflow_error Invalid numeric cast (underflow).
*/
template< typename To, typename From >
inline EnableIf_< IsCriticalFloatFloatConversion<To,From>, To >
   numeric_cast_backend( From from )
{
   if( from > From( std::numeric_limits<To>::max() ) )
      BLAZE_THROW_OVERFLOW_ERROR( "Invalid numeric cast (overflow)" );

   if( from < From( std::numeric_limits<To>::min() ) )
      BLAZE_THROW_UNDERFLOW_ERROR( "Invalid numeric cast (underflow)" );

   BLAZE_INTERNAL_ASSERT( from == From( To( from ) ), "Numeric cast failed" );

   return from;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked conversion of values of numeric type.
// \ingroup util
//
// \param from The numeric value to be converted.
// \return The converted value.
// \exception std::overflow_error Invalid numeric cast (overflow).
// \exception std::underflow_error Invalid numeric cast (underflow).
//
// This function converts the given numeric value \a from to the specified type \a To. In case
// a loss of range is detected, either a \a std::underflow_error or \a std::overflow_error
// exception is thrown.
//
// Examples:

   \code
   // Triggers a std::overflow_error exception
   try {
      const int a( std::numeric_limits<int>::max() );
      const short b( blaze::numeric_cast<short>( a ) );
   }
   catch( std::overflow_error& ) {}

   // Triggers a std::underflow_error exception
   try {
      const int a( -1 );
      const unsigned int b( blaze::numeric_cast<unsigned int>( a ) );
   }
   catch( std::underflow_error& ) {}
   \endcode
*/
template< typename To      // The target type
        , typename From >  // The source type
inline To numeric_cast( From from )
{
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( To   );
   BLAZE_CONSTRAINT_MUST_BE_NUMERIC_TYPE( From );

   return numeric_cast_backend<To>( from );
}
//*************************************************************************************************

} // namespace blaze

#endif
