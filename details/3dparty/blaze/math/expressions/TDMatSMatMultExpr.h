//=================================================================================================
/*!
//  \file blaze/math/expressions/TDMatSMatMultExpr.h
//  \brief Header file for the transpose dense matrix/sparse matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TDMATSMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDMATSMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "../../math/Aliases.h"
#include "../../math/constraints/ColumnMajorMatrix.h"
#include "../../math/constraints/DenseMatrix.h"
#include "../../math/constraints/MatMatMultExpr.h"
#include "../../math/constraints/RequiresEvaluation.h"
#include "../../math/constraints/RowMajorMatrix.h"
#include "../../math/constraints/SparseMatrix.h"
#include "../../math/constraints/StorageOrder.h"
#include "../../math/constraints/Symmetric.h"
#include "../../math/Exception.h"
#include "../../math/expressions/Computation.h"
#include "../../math/expressions/DenseMatrix.h"
#include "../../math/expressions/Forward.h"
#include "../../math/expressions/MatMatMultExpr.h"
#include "../../math/functors/DeclDiag.h"
#include "../../math/functors/DeclHerm.h"
#include "../../math/functors/DeclLow.h"
#include "../../math/functors/DeclSym.h"
#include "../../math/functors/DeclUpp.h"
#include "../../math/functors/Noop.h"
#include "../../math/shims/Conjugate.h"
#include "../../math/shims/IsDefault.h"
#include "../../math/shims/Reset.h"
#include "../../math/shims/Serial.h"
#include "../../math/traits/MultTrait.h"
#include "../../math/typetraits/IsAligned.h"
#include "../../math/typetraits/IsColumnMajorMatrix.h"
#include "../../math/typetraits/IsComputation.h"
#include "../../math/typetraits/IsDiagonal.h"
#include "../../math/typetraits/IsExpression.h"
#include "../../math/typetraits/IsHermitian.h"
#include "../../math/typetraits/IsLower.h"
#include "../../math/typetraits/IsResizable.h"
#include "../../math/typetraits/IsStrictlyLower.h"
#include "../../math/typetraits/IsStrictlyUpper.h"
#include "../../math/typetraits/IsSymmetric.h"
#include "../../math/typetraits/IsTriangular.h"
#include "../../math/typetraits/IsUniLower.h"
#include "../../math/typetraits/IsUniUpper.h"
#include "../../math/typetraits/IsUpper.h"
#include "../../math/typetraits/RequiresEvaluation.h"
#include "../../math/typetraits/Size.h"
#include "../../math/views/Check.h"
#include "../../system/Optimizations.h"
#include "../../system/Thresholds.h"
#include "../../util/algorithms/Max.h"
#include "../../util/algorithms/Min.h"
#include "../../util/Assert.h"
#include "../../util/DisableIf.h"
#include "../../util/EnableIf.h"
#include "../../util/FunctionTrace.h"
#include "../../util/IntegralConstant.h"
#include "../../util/mpl/If.h"
#include "../../util/Types.h"
#include "../../util/typetraits/IsBuiltin.h"


namespace blaze {

//=================================================================================================
//
//  CLASS TDMATSMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose dense matrix-sparse matrix multiplications.
// \ingroup dense_matrix_expression
//
// The TDMatSMatMultExpr class represents the compile time expression for multiplications between
// a column-major dense matrix and a row-major sparse matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
class TDMatSMatMultExpr
   : public MatMatMultExpr< DenseMatrix< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF>, true > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<MT1>;     //!< Result type of the left-hand side dense matrix expression.
   using RT2 = ResultType_t<MT2>;     //!< Result type of the right-hand side sparse matrix expression.
   using ET1 = ElementType_t<RT1>;    //!< Element type of the left-hand side dense matrix expression.
   using ET2 = ElementType_t<RT2>;    //!< Element type of the right-hand side sparse matrix expression.
   using CT1 = CompositeType_t<MT1>;  //!< Composite type of the left-hand side dense matrix expression.
   using CT2 = CompositeType_t<MT2>;  //!< Composite type of the right-hand side sparse matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense matrix expression.
   static constexpr bool evaluateLeft = ( IsComputation_v<MT1> || RequiresEvaluation_v<MT1> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side sparse matrix expression.
   static constexpr bool evaluateRight = ( IsComputation_v<MT2> || RequiresEvaluation_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   static constexpr bool SYM  = ( SF && !( HF || LF || UF )    );  //!< Flag for symmetric matrices.
   static constexpr bool HERM = ( HF && !( LF || UF )          );  //!< Flag for Hermitian matrices.
   static constexpr bool LOW  = ( LF || ( ( SF || HF ) && UF ) );  //!< Flag for lower matrices.
   static constexpr bool UPP  = ( UF || ( ( SF || HF ) && LF ) );  //!< Flag for upper matrices.
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the optimal evaluation strategy.
       In case the right-hand side matrix operands of type \a T3 is symmetric, the variable is
       set to 1 and an optimized evaluation strategy is selected. Otherwise the variable is set
       to 0 and the default strategy is chosen. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool CanExploitSymmetry_v = IsSymmetric_v<T3>;
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either of the two matrix operands requires an intermediate evaluation, the variable
       will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool IsEvaluationRequired_v =
      ( ( evaluateLeft || evaluateRight ) && !CanExploitSymmetry_v<T1,T2,T3> );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case no SMP assignment is required and the element type of the target matrix has a
       fixed size (i.e. is not resizable), the variable will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseOptimizedKernel_v =
      ( useOptimizedKernels &&
        !IsDiagonal_v<T2> &&
        !IsResizable_v< ElementType_t<T1> > &&
        !( IsColumnMajorMatrix_v<T1> && IsResizable_v<ET2> ) );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case no SMP assignment is required and the element type of the target matrix is
       resizable, the variable will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseDefaultKernel_v = !UseOptimizedKernel_v<T1,T2,T3>;
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Type of the functor for forwarding an expression to another assign kernel.
   /*! In case a temporary matrix needs to be created, this functor is used to forward the
       resulting expression to another assign kernel. */
   using ForwardFunctor = If_t< HERM
                              , DeclHerm
                              , If_t< SYM
                                    , DeclSym
                                    , If_t< LOW
                                          , If_t< UPP
                                                , DeclDiag
                                                , DeclLow >
                                          , If_t< UPP
                                                , DeclUpp
                                                , Noop > > > >;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this TDMatSMatMultExpr instance.
   using This = TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF>;

   using ResultType    = MultTrait_t<RT1,RT2>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT1>, const MT1, const MT1& >;

   //! Composite type of the right-hand side sparse matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;

   //! Type for the assignment of the left-hand side dense matrix operand.
   using LT = If_t< evaluateLeft, const RT1, CT1 >;

   //! Type for the assignment of the right-hand side sparse matrix operand.
   using RT = If_t< evaluateRight, const RT2, CT2 >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable =
      ( !evaluateLeft && MT1::smpAssignable && !evaluateRight && MT2::smpAssignable );
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TDMatSMatMultExpr class.
   //
   // \param lhs The left-hand side dense matrix operand of the multiplication expression.
   // \param rhs The right-hand side sparse matrix operand of the multiplication expression.
   */
   explicit inline TDMatSMatMultExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side dense matrix of the multiplication expression
      , rhs_( rhs )  // Right-hand side sparse matrix of the multiplication expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.columns() == rhs.rows(), "Invalid matrix sizes" );
   }
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const {
      BLAZE_INTERNAL_ASSERT( i < lhs_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < rhs_.columns(), "Invalid column access index" );

      if( IsDiagonal_v<MT1> ) {
         return lhs_(i,i) * rhs_(i,j);
      }
      else if( IsDiagonal_v<MT2> ) {
         return lhs_(i,j) * rhs_(j,j);
      }
      else if( IsTriangular_v<MT1> || IsTriangular_v<MT2> ) {
         const size_t begin( ( IsUpper_v<MT1> )
                             ?( ( IsLower_v<MT2> )
                                ?( max( ( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                                      , ( IsStrictlyLower_v<MT2> ? j+1UL : j ) ) )
                                :( IsStrictlyUpper_v<MT1> ? i+1UL : i ) )
                             :( ( IsLower_v<MT2> )
                                ?( IsStrictlyLower_v<MT2> ? j+1UL : j )
                                :( 0UL ) ) );
         const size_t end( ( IsLower_v<MT1> )
                           ?( ( IsUpper_v<MT2> )
                              ?( min( ( IsStrictlyLower_v<MT1> ? i : i+1UL )
                                    , ( IsStrictlyUpper_v<MT2> ? j : j+1UL ) ) )
                              :( IsStrictlyLower_v<MT1> ? i : i+1UL ) )
                           :( ( IsUpper_v<MT2> )
                              ?( IsStrictlyUpper_v<MT2> ? j : j+1UL )
                              :( lhs_.columns() ) ) );

         if( begin >= end ) return ElementType();

         const size_t n( end - begin );

         return subvector( row( lhs_, i, unchecked ), begin, n, unchecked ) *
                subvector( column( rhs_, j, unchecked ), begin, n, unchecked );
      }
      else {
         return row( lhs_, i, unchecked ) * column( rhs_, j, unchecked );
      }
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid matrix access index.
   */
   inline ReturnType at( size_t i, size_t j ) const {
      if( i >= lhs_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= rhs_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return lhs_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return rhs_.columns();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side transpose dense matrix operand.
   //
   // \return The left-hand side transpose dense matrix operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side sparse matrix operand.
   //
   // \return The right-hand side sparse matrix operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return rhs_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const {
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return lhs_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return ( rows() * columns() >= SMP_TDMATSMATMULT_THRESHOLD ) && !IsDiagonal_v<MT1>;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side sparse matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-sparse matrix multiplication to a dense matrix
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      assign( DenseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TDMatSMatMultExpr::selectAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to dense matrices (kernel selection)*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for an assignment of a transpose dense matrix-sparse matrix
   //        multiplication to a dense matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      if( C.rows() * C.columns() < TDMATSMATMULT_THRESHOLD )
         selectSmallAssignKernel( C, A, B );
      else
         selectLargeAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default assignment of a transpose dense matrix-sparse matrix
   // multiplication expression to a dense matrix. This assign function is used in case the
   // element type of the target matrix is resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectDefaultAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      using ConstIterator = ConstIterator_t<MT5>;

      reset( C );

      for( size_t j=0UL; j<B.rows(); ++j )
      {
         ConstIterator element( B.begin(j) );
         const ConstIterator end( B.end(j) );

         if( IsDiagonal_v<MT4> )
         {
            for( ; element!=end; ++element ) {
               C(j,element->index()) = A(j,j) * element->value();
            }
         }
         else
         {
            for( ; element!=end; ++element )
            {
               const size_t j1( element->index() );

               const size_t ibegin( ( IsLower_v<MT4> )
                                    ?( ( IsStrictlyLower_v<MT4> )
                                       ?( SYM || HERM || LOW ? max(j1,j+1UL) : j+1UL )
                                       :( SYM || HERM || LOW ? max(j1,j) : j ) )
                                    :( SYM || HERM || LOW ? j1 : 0UL ) );
               const size_t iend( ( IsUpper_v<MT4> )
                                  ?( ( IsStrictlyUpper_v<MT4> )
                                     ?( UPP ? min(j1+1UL,j) : j )
                                     :( UPP ? min(j1,j)+1UL : j+1UL ) )
                                  :( UPP ? j1+1UL : A.rows() ) );

               if( ( SYM || HERM || LOW || UPP ) && ( ibegin >= iend ) ) continue;
               BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

               for( size_t i=ibegin; i<iend; ++i ) {
                  if( isDefault( C(i,j1) ) )
                     C(i,j1) = A(i,j) * element->value();
                  else
                     C(i,j1) += A(i,j) * element->value();
               }
            }
         }
      }

      if( SYM || HERM ) {
         for( size_t j=1UL; j<B.columns(); ++j ) {
            for( size_t i=0UL; i<j; ++i ) {
               C(i,j) = HERM ? conj( C(j,i) ) : C(j,i);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices (small matrices)***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a small transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseDefaultKernel_v<MT3,MT4,MT5> >
      selectSmallAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      selectDefaultAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to dense matrices (small matrices)*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized assignment of a small transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a dense matrix. This assign function is used in
   // case the element type of the target matrix is not resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseOptimizedKernel_v<MT3,MT4,MT5> >
      selectSmallAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      using ConstIterator = ConstIterator_t<MT5>;

      reset( C );

      for( size_t j=0UL; j<B.rows(); ++j )
      {
         ConstIterator element( B.begin(j) );
         const ConstIterator end( B.end(j) );

         const size_t nonzeros( B.nonZeros(j) );
         const size_t kpos( nonzeros & size_t(-4) );
         BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kpos, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );
            ++element;
            const size_t j2( element->index() );
            const ET2    v2( element->value() );
            ++element;
            const size_t j3( element->index() );
            const ET2    v3( element->value() );
            ++element;
            const size_t j4( element->index() );
            const ET2    v4( element->value() );
            ++element;

            BLAZE_INTERNAL_ASSERT( j1 < j2 && j2 < j3 && j3 < j4, "Invalid sparse matrix index detected" );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( SYM || HERM || LOW ? max(j1,j+1UL) : j+1UL )
                                    :( SYM || HERM || LOW ? max(j1,j) : j ) )
                                 :( SYM || HERM || LOW ? j1 : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j4+1UL,j) : j )
                                  :( UPP ? min(j4,j)+1UL : j+1UL ) )
                               :( UPP ? j4+1UL : A.rows() ) );

            if( ( SYM || HERM || LOW || UPP ) && ( ibegin >= iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( ibegin + inum - ( inum % 4UL ) ) == ipos, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=4UL ) {
               C(i    ,j1) += A(i    ,j) * v1;
               C(i+1UL,j1) += A(i+1UL,j) * v1;
               C(i+2UL,j1) += A(i+2UL,j) * v1;
               C(i+3UL,j1) += A(i+3UL,j) * v1;
               C(i    ,j2) += A(i    ,j) * v2;
               C(i+1UL,j2) += A(i+1UL,j) * v2;
               C(i+2UL,j2) += A(i+2UL,j) * v2;
               C(i+3UL,j2) += A(i+3UL,j) * v2;
               C(i    ,j3) += A(i    ,j) * v3;
               C(i+1UL,j3) += A(i+1UL,j) * v3;
               C(i+2UL,j3) += A(i+2UL,j) * v3;
               C(i+3UL,j3) += A(i+3UL,j) * v3;
               C(i    ,j4) += A(i    ,j) * v4;
               C(i+1UL,j4) += A(i+1UL,j) * v4;
               C(i+2UL,j4) += A(i+2UL,j) * v4;
               C(i+3UL,j4) += A(i+3UL,j) * v4;
            }
            for( ; i<iend; ++i ) {
               C(i,j1) += A(i,j) * v1;
               C(i,j2) += A(i,j) * v2;
               C(i,j3) += A(i,j) * v3;
               C(i,j4) += A(i,j) * v4;
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( SYM || HERM || LOW ? max(j1,j+1UL) : j+1UL )
                                    :( SYM || HERM || LOW ? max(j1,j) : j ) )
                                 :( SYM || HERM || LOW ? j1 : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j1+1UL,j) : j )
                                  :( UPP ? min(j1,j)+1UL : j+1UL ) )
                               :( UPP ? j1+1UL : A.rows() ) );

            if( ( SYM || HERM || LOW || UPP ) && ( ibegin >= iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( ibegin + inum - ( inum % 4UL ) ) == ipos, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=4UL ) {
               C(i    ,j1) += A(i    ,j) * v1;
               C(i+1UL,j1) += A(i+1UL,j) * v1;
               C(i+2UL,j1) += A(i+2UL,j) * v1;
               C(i+3UL,j1) += A(i+3UL,j) * v1;
            }
            for( ; i<iend; ++i ) {
               C(i,j1) += A(i,j) * v1;
            }
         }
      }

      if( SYM || HERM ) {
         for( size_t j=1UL; j<B.columns(); ++j ) {
            for( size_t i=0UL; i<j; ++i ) {
               C(i,j) = HERM ? conj( C(j,i) ) : C(j,i);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices (large matrices)***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a large transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseDefaultKernel_v<MT3,MT4,MT5> >
      selectLargeAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      selectDefaultAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to dense matrices (large matrices)*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a large transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a dense matrix. This assign function is used in
   // case the element type of the target matrix is not resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseOptimizedKernel_v<MT3,MT4,MT5> >
      selectLargeAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OppositeType_t<MT5> );

      const ForwardFunctor fwd;

      const OppositeType_t<MT5> tmp( serial( B ) );
      assign( C, fwd( A * tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-sparse matrix multiplication to a sparse matrix
   //        (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      assign( SparseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      const TmpType tmp( serial( rhs ) );
      assign( ~lhs, fwd( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring assignment********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring assignment of a transpose dense
   // matrix-sparse matrix multiplication expression. Due to the explicit application of the
   // SFINAE principle this function can only be selected by the compiler in case the symmetry
   // of either of the two matrix operands can be exploited.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      assign( Matrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      assign( ~lhs, fwd( rhs.lhs_ * trans( rhs.rhs_ ) ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose dense matrix-sparse matrix multiplication to
   //        a dense matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      addAssign( DenseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TDMatSMatMultExpr::selectAddAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices (kernel selection)************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for an addition assignment of a transpose dense matrix-sparse
   //        matrix multiplication to a dense matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      if( C.rows() * C.columns() < TDMATSMATMULT_THRESHOLD )
         selectSmallAddAssignKernel( C, A, B );
      else
         selectLargeAddAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices***********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default addition assignment of a transpose dense matrix-sparse
   // matrix multiplication expression to a dense matrix. This assign function is use in case the
   // element type of the target matrix is resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectDefaultAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      using ConstIterator = ConstIterator_t<MT5>;

      size_t i( 0UL );

      for( size_t j=0UL; j<B.rows(); ++j )
      {
         ConstIterator element( B.begin(j) );
         const ConstIterator end( B.end(j) );

         if( IsDiagonal_v<MT4> )
         {
            for( ; element!=end; ++element ) {
               C(j,element->index()) += A(j,j) * element->value();
            }
         }
         else
         {
            for( ; element!=end; ++element )
            {
               const size_t j1( element->index() );

               const size_t ibegin( ( IsLower_v<MT4> )
                                    ?( ( IsStrictlyLower_v<MT4> )
                                       ?( LOW ? max(j1,j+1UL) : j+1UL )
                                       :( LOW ? max(j1,j) : j ) )
                                    :( LOW ? j1 : 0UL ) );
               const size_t iend( ( IsUpper_v<MT4> )
                                  ?( ( IsStrictlyUpper_v<MT4> )
                                     ?( UPP ? min(j1+1UL,j) : j )
                                     :( UPP ? min(j1,j)+1UL : j+1UL ) )
                                  :( UPP ? j1+1UL : A.rows() ) );

               if( ( LOW || UPP ) && ( ibegin >= iend ) ) continue;
               BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

               const size_t inum( iend - ibegin );
               const size_t ipos( ibegin + ( inum & size_t(-4) ) );
               BLAZE_INTERNAL_ASSERT( ( ibegin + inum - ( inum % 4UL ) ) == ipos, "Invalid end calculation" );

               for( i=ibegin; i<ipos; i+=4UL ) {
                  C(i    ,j1) += A(i    ,j) * element->value();
                  C(i+1UL,j1) += A(i+1UL,j) * element->value();
                  C(i+2UL,j1) += A(i+2UL,j) * element->value();
                  C(i+3UL,j1) += A(i+3UL,j) * element->value();
               }
               for( ; i<iend; ++i ) {
                  C(i,j1) += A(i,j) * element->value();
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (small matrices)******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a small transpose dense matrix-sparse matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseDefaultKernel_v<MT3,MT4,MT5> >
      selectSmallAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      selectDefaultAddAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized addition assignment to dense matrices (small matrices)****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized addition assignment of a small transpose dense matrix-sparse matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a dense matrix. This assign function is
   // used in case the element type of the target matrix is not resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseOptimizedKernel_v<MT3,MT4,MT5> >
      selectSmallAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      using ConstIterator = ConstIterator_t<MT5>;

      for( size_t j=0UL; j<B.rows(); ++j )
      {
         ConstIterator element( B.begin(j) );
         const ConstIterator end( B.end(j) );

         const size_t nonzeros( B.nonZeros(j) );
         const size_t kpos( nonzeros & size_t(-4) );
         BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kpos, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );
            ++element;
            const size_t j2( element->index() );
            const ET2    v2( element->value() );
            ++element;
            const size_t j3( element->index() );
            const ET2    v3( element->value() );
            ++element;
            const size_t j4( element->index() );
            const ET2    v4( element->value() );
            ++element;

            BLAZE_INTERNAL_ASSERT( j1 < j2 && j2 < j3 && j3 < j4, "Invalid sparse matrix index detected" );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( LOW ? max(j1,j+1UL) : j+1UL )
                                    :( LOW ? max(j1,j) : j ) )
                                 :( LOW ? j1 : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j4+1UL,j) : j )
                                  :( UPP ? min(j4,j)+1UL : j+1UL ) )
                               :( UPP ? j4+1UL : A.rows() ) );

            if( ( LOW || UPP ) && ( ibegin >= iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( ibegin + inum - ( inum % 4UL ) ) == ipos, "Invalid end calculation" );

            size_t i( ibegin );

            for( i=ibegin; i<ipos; i+=4UL ) {
               C(i    ,j1) += A(i    ,j) * v1;
               C(i+1UL,j1) += A(i+1UL,j) * v1;
               C(i+2UL,j1) += A(i+2UL,j) * v1;
               C(i+3UL,j1) += A(i+3UL,j) * v1;
               C(i    ,j2) += A(i    ,j) * v2;
               C(i+1UL,j2) += A(i+1UL,j) * v2;
               C(i+2UL,j2) += A(i+2UL,j) * v2;
               C(i+3UL,j2) += A(i+3UL,j) * v2;
               C(i    ,j3) += A(i    ,j) * v3;
               C(i+1UL,j3) += A(i+1UL,j) * v3;
               C(i+2UL,j3) += A(i+2UL,j) * v3;
               C(i+3UL,j3) += A(i+3UL,j) * v3;
               C(i    ,j4) += A(i    ,j) * v4;
               C(i+1UL,j4) += A(i+1UL,j) * v4;
               C(i+2UL,j4) += A(i+2UL,j) * v4;
               C(i+3UL,j4) += A(i+3UL,j) * v4;
            }
            for( ; i<iend; ++i ) {
               C(i,j1) += A(i,j) * v1;
               C(i,j2) += A(i,j) * v2;
               C(i,j3) += A(i,j) * v3;
               C(i,j4) += A(i,j) * v4;
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( LOW ? max(j1,j+1UL) : j+1UL )
                                    :( LOW ? max(j1,j) : j ) )
                                 :( LOW ? j1 : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j1+1UL,j) : j )
                                  :( UPP ? min(j1,j)+1UL : j+1UL ) )
                               :( UPP ? j1+1UL : A.rows() ) );

            if( ( LOW || UPP ) && ( ibegin >= iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( ibegin + inum - ( inum % 4UL ) ) == ipos, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=4UL ) {
               C(i    ,j1) += A(i    ,j) * v1;
               C(i+1UL,j1) += A(i+1UL,j) * v1;
               C(i+2UL,j1) += A(i+2UL,j) * v1;
               C(i+3UL,j1) += A(i+3UL,j) * v1;
            }
            for( ; i<iend; ++i ) {
               C(i,j1) += A(i,j) * v1;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (large matrices)******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a large transpose dense matrix-sparse matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseDefaultKernel_v<MT3,MT4,MT5> >
      selectLargeAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      selectDefaultAddAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized addition assignment to dense matrices (large matrices)****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a large transpose dense matrix-sparse matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a dense matrix. This assign function is
   // used in case the element type of the target matrix is not resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseOptimizedKernel_v<MT3,MT4,MT5> >
      selectLargeAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OppositeType_t<MT5> );

      const ForwardFunctor fwd;

      const OppositeType_t<MT5> tmp( serial( B ) );
      addAssign( C, fwd( A * tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring addition assignment***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring addition assignment of a transpose dense matrix-sparse matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring addition assignment of a transpose
   // dense matrix-sparse matrix multiplication expression. Due to the explicit application of the
   // SFINAE principle this function can only be selected by the compiler in case the symmetry of
   // either of the two matrix operands can be exploited.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      addAssign( Matrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      addAssign( ~lhs, fwd( rhs.lhs_ * trans( rhs.rhs_ ) ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose dense matrix-sparse matrix multiplication to a
   //        dense matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      subAssign( DenseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      TDMatSMatMultExpr::selectSubAssignKernel( ~lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices (kernel selection)*********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for a subtraction assignment of a transpose dense matrix-
   //        sparse matrix multiplication to a dense matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      if( C.rows() * C.columns() < TDMATSMATMULT_THRESHOLD )
         selectSmallSubAssignKernel( C, A, B );
      else
         selectLargeSubAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose dense matrix-sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default subtraction assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a dense matrix. This assign function is used
   // in case the element type of the target matrix is resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectDefaultSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      using ConstIterator = ConstIterator_t<MT5>;

      size_t i( 0UL );

      for( size_t j=0UL; j<B.rows(); ++j )
      {
         ConstIterator element( B.begin(j) );
         const ConstIterator end( B.end(j) );

         if( IsDiagonal_v<MT4> )
         {
            for( ; element!=end; ++element ) {
               C(j,element->index()) -= A(j,j) * element->value();
            }
         }
         else
         {
            for( ; element!=end; ++element )
            {
               const size_t j1( element->index() );

               const size_t ibegin( ( IsLower_v<MT4> )
                                    ?( ( IsStrictlyLower_v<MT4> )
                                       ?( LOW ? max(j1,j+1UL) : j+1UL )
                                       :( LOW ? max(j1,j) : j ) )
                                    :( LOW ? j1 : 0UL ) );
               const size_t iend( ( IsUpper_v<MT4> )
                                  ?( ( IsStrictlyUpper_v<MT4> )
                                     ?( UPP ? min(j1+1UL,j) : j )
                                     :( UPP ? min(j1,j)+1UL : j+1UL ) )
                                  :( UPP ? j1+1UL : A.rows() ) );

               if( ( LOW || UPP ) && ( ibegin >= iend ) ) continue;
               BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

               const size_t inum( iend - ibegin );
               const size_t ipos( ibegin + ( inum & size_t(-4) ) );
               BLAZE_INTERNAL_ASSERT( ( ibegin + inum - ( inum % 4UL ) ) == ipos, "Invalid end calculation" );

               for( i=ibegin; i<ipos; i+=4UL ) {
                  C(i    ,j1) -= A(i    ,j) * element->value();
                  C(i+1UL,j1) -= A(i+1UL,j) * element->value();
                  C(i+2UL,j1) -= A(i+2UL,j) * element->value();
                  C(i+3UL,j1) -= A(i+3UL,j) * element->value();
               }
               for( ; i<iend; ++i ) {
                  C(i,j1) -= A(i,j) * element->value();
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (small matrices)***************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a small transpose dense matrix-sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a
   // transpose dense matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseDefaultKernel_v<MT3,MT4,MT5> >
      selectSmallSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      selectDefaultSubAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized subtraction assignment to dense matrices (small matrices)*************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized subtraction assignment of a small transpose dense matrix-sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a dense matrix. This assign function
   // is used in case the element type of the target matrix is not resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseOptimizedKernel_v<MT3,MT4,MT5> >
      selectSmallSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      using ConstIterator = ConstIterator_t<MT5>;

      for( size_t j=0UL; j<B.rows(); ++j )
      {
         ConstIterator element( B.begin(j) );
         const ConstIterator end( B.end(j) );

         const size_t nonzeros( B.nonZeros(j) );
         const size_t kpos( nonzeros & size_t(-4) );
         BLAZE_INTERNAL_ASSERT( ( nonzeros - ( nonzeros % 4UL ) ) == kpos, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );
            ++element;
            const size_t j2( element->index() );
            const ET2    v2( element->value() );
            ++element;
            const size_t j3( element->index() );
            const ET2    v3( element->value() );
            ++element;
            const size_t j4( element->index() );
            const ET2    v4( element->value() );
            ++element;

            BLAZE_INTERNAL_ASSERT( j1 < j2 && j2 < j3 && j3 < j4, "Invalid sparse matrix index detected" );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( LOW ? max(j1,j+1UL) : j+1UL )
                                    :( LOW ? max(j1,j) : j ) )
                                 :( LOW ? j1 : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j4+1UL,j) : j )
                                  :( UPP ? min(j4,j)+1UL : j+1UL ) )
                               :( UPP ? j4+1UL : A.rows() ) );

            if( ( LOW || UPP ) && ( ibegin >= iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( ibegin + inum - ( inum % 4UL ) ) == ipos, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=4UL ) {
               C(i    ,j1) -= A(i    ,j) * v1;
               C(i+1UL,j1) -= A(i+1UL,j) * v1;
               C(i+2UL,j1) -= A(i+2UL,j) * v1;
               C(i+3UL,j1) -= A(i+3UL,j) * v1;
               C(i    ,j2) -= A(i    ,j) * v2;
               C(i+1UL,j2) -= A(i+1UL,j) * v2;
               C(i+2UL,j2) -= A(i+2UL,j) * v2;
               C(i+3UL,j2) -= A(i+3UL,j) * v2;
               C(i    ,j3) -= A(i    ,j) * v3;
               C(i+1UL,j3) -= A(i+1UL,j) * v3;
               C(i+2UL,j3) -= A(i+2UL,j) * v3;
               C(i+3UL,j3) -= A(i+3UL,j) * v3;
               C(i    ,j4) -= A(i    ,j) * v4;
               C(i+1UL,j4) -= A(i+1UL,j) * v4;
               C(i+2UL,j4) -= A(i+2UL,j) * v4;
               C(i+3UL,j4) -= A(i+3UL,j) * v4;
            }
            for( ; i<iend; ++i ) {
               C(i,j1) -= A(i,j) * v1;
               C(i,j2) -= A(i,j) * v2;
               C(i,j3) -= A(i,j) * v3;
               C(i,j4) -= A(i,j) * v4;
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( LOW ? max(j1,j+1UL) : j+1UL )
                                    :( LOW ? max(j1,j) : j ) )
                                 :( LOW ? j1 : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j1+1UL,j) : j )
                                  :( UPP ? min(j1,j)+1UL : j+1UL ) )
                               :( UPP ? j1+1UL : A.rows() ) );

            if( ( LOW || UPP ) && ( ibegin >= iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + ( inum & size_t(-4) ) );
            BLAZE_INTERNAL_ASSERT( ( ibegin + inum - ( inum % 4UL ) ) == ipos, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=4UL ) {
               C(i    ,j1) -= A(i    ,j) * v1;
               C(i+1UL,j1) -= A(i+1UL,j) * v1;
               C(i+2UL,j1) -= A(i+2UL,j) * v1;
               C(i+3UL,j1) -= A(i+3UL,j) * v1;
            }
            for( ; i<iend; ++i ) {
               C(i,j1) -= A(i,j) * v1;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (large matrices)***************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a large transpose dense matrix-sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a
   // transpose dense matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseDefaultKernel_v<MT3,MT4,MT5> >
      selectLargeSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      selectDefaultSubAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized subtraction assignment to dense matrices (large matrices)*************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a large transpose dense matrix-sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a dense matrix. This assign function
   // is used in case the element type of the target matrix is not resizable.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline EnableIf_t< UseOptimizedKernel_v<MT3,MT4,MT5> >
      selectLargeSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OppositeType_t<MT5> );

      const ForwardFunctor fwd;

      const OppositeType_t<MT5> tmp( serial( B ) );
      subAssign( C, fwd( A * tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring subtraction assignment********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring subtraction assignment of a transpose dense matrix-sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring subtraction assignment of a
   // transpose dense matrix-sparse matrix multiplication expression. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler
   // in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      subAssign( Matrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      subAssign( ~lhs, fwd( rhs.lhs_ * trans( rhs.rhs_ ) ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a transpose dense matrix-sparse matrix multiplication
   //        to a dense matrix (\f$ C\circ=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      schurAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to sparse matrices*************************************************
   // No special implementation for the Schur product assignment to sparse matrices.
   //**********************************************************************************************

   //**Multiplication assignment to dense matrices*************************************************
   // No special implementation for the multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**Multiplication assignment to sparse matrices************************************************
   // No special implementation for the multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-sparse matrix multiplication to a dense matrix
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // sparse matrix multiplication expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle this function can only be selected by the compiler in case either
   // of the two matrix operands requires an intermediate evaluation and no symmetry can be
   // exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
      smpAssign( DenseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      smpAssign( ~lhs, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose dense matrix-sparse matrix multiplication to a sparse
   //        matrix (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a sparse matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler
   // in case either of the two matrix operands requires an intermediate evaluation and no
   // symmetry can be exploited.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
      smpAssign( SparseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      const TmpType tmp( rhs );
      smpAssign( ~lhs, fwd( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring SMP assignment****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP assignment of a transpose dense matrix-sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP assignment of a transpose
   // dense matrix-sparse matrix multiplication expression. Due to the explicit application
   // of the SFINAE principle this function can only be selected by the compiler in case the
   // symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      smpAssign( Matrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      smpAssign( ~lhs, fwd( rhs.lhs_ * trans( rhs.rhs_ ) ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose dense matrix-sparse matrix multiplication to
   //        a dense matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-sparse matrix multiplication expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler in
   // case either of the two matrix operands requires an intermediate evaluation and no symmetry
   // can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
      smpAddAssign( DenseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      smpAddAssign( ~lhs, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring SMP addition assignment*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP addition assignment of a transpose dense matrix-sparse matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP addition assignment of
   // a transpose dense matrix-sparse matrix multiplication expression. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler
   // in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      smpAddAssign( Matrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      smpAddAssign( ~lhs, fwd( rhs.lhs_ * trans( rhs.rhs_ ) ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose dense matrix-sparse matrix multiplication
   //        to a dense matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a dense matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler in
   // case either of the two matrix operands requires an intermediate evaluation and no symmetry
   // can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
      smpSubAssign( DenseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (~lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (~lhs).columns()  , "Invalid number of columns" );

      smpSubAssign( ~lhs, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring SMP subtraction assignment****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP subtraction assignment of a transpose dense matrix-sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP subtraction assignment of
   // a transpose dense matrix-sparse matrix multiplication expression. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler
   // in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
      smpSubAssign( Matrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      smpSubAssign( ~lhs, fwd( rhs.lhs_ * trans( rhs.rhs_ ) ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a transpose dense matrix-sparse matrix multiplication
   //        to a dense matrix (\f$ C\circ=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a transpose
   // dense matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline void smpSchurAssign( DenseMatrix<MT,SO>& lhs, const TDMatSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (~lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (~lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      smpSchurAssign( ~lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP Schur product assignment to sparse matrices*********************************************
   // No special implementation for the SMP Schur product assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense matrices*********************************************
   // No special implementation for the SMP multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse matrices********************************************
   // No special implementation for the SMP multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATMATMULTEXPR( MT1, MT2 );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a column-major dense matrix and a
//        row-major sparse matrix (\f$ A=B*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the multiplication.
// \param rhs The right-hand side sparse matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of a column-major dense matrix and a row-major
// sparse matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,rowMajor> B;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline decltype(auto)
   operator*( const DenseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (~lhs).columns() != (~rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }


   using ReturnType = const TDMatSMatMultExpr<MT1,MT2,false,false,false,false>;
   return ReturnType( ~lhs, ~rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-symmetric matrix multiplication expression as symmetric.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid symmetric matrix specification.
//
// The \a declsym function declares the given non-symmetric matrix multiplication expression
// \a dm as symmetric. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument exception
// is thrown.\n
// The following example demonstrates the use of the \a declsym function:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,rowMajor> B;
   // ... Resizing and initialization
   C = declsym( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) declsym( const TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid symmetric matrix specification" );
   }

   using ReturnType = const TDMatSMatMultExpr<MT1,MT2,true,HF,LF,UF>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-Hermitian matrix multiplication expression as Hermitian.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid Hermitian matrix specification.
//
// The \a declherm function declares the given non-Hermitian matrix multiplication expression
// \a dm as Hermitian. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument exception
// is thrown.\n
// The following example demonstrates the use of the \a declherm function:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,rowMajor> B;
   // ... Resizing and initialization
   C = declherm( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) declherm( const TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid Hermitian matrix specification" );
   }

   using ReturnType = const TDMatSMatMultExpr<MT1,MT2,SF,true,LF,UF>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-lower matrix multiplication expression as lower.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid lower matrix specification.
//
// The \a decllow function declares the given non-lower matrix multiplication expression
// \a dm as lower. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument
// exception is thrown.\n
// The following example demonstrates the use of the \a decllow function:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,rowMajor> B;
   // ... Resizing and initialization
   C = decllow( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) decllow( const TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid lower matrix specification" );
   }

   using ReturnType = const TDMatSMatMultExpr<MT1,MT2,SF,HF,true,UF>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-upper matrix multiplication expression as upper.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid upper matrix specification.
//
// The \a declupp function declares the given non-upper matrix multiplication expression
// \a dm as upper. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument
// exception is thrown.\n
// The following example demonstrates the use of the \a declupp function:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,rowMajor> B;
   // ... Resizing and initialization
   C = declupp( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) declupp( const TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid upper matrix specification" );
   }

   using ReturnType = const TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,true>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-diagonal matrix multiplication expression as diagonal.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid diagonal matrix specification.
//
// The \a decldiag function declares the given non-diagonal matrix multiplication expression
// \a dm as diagonal. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument exception
// is thrown.\n
// The following example demonstrates the use of the \a decldiag function:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,rowMajor> B;
   // ... Resizing and initialization
   C = decldiag( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) decldiag( const TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid diagonal matrix specification" );
   }

   using ReturnType = const TDMatSMatMultExpr<MT1,MT2,SF,HF,true,true>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct Size< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF>, 0UL >
   : public Size<MT1,0UL>
{};

template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct Size< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF>, 1UL >
   : public Size<MT2,1UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct IsAligned< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public IsAligned<MT1>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct IsSymmetric< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public BoolConstant< SF ||
                          ( HF && IsBuiltin_v< ElementType_t< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF> > > ) ||
                          ( LF && UF ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISHERMITIAN SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool LF, bool UF >
struct IsHermitian< TDMatSMatMultExpr<MT1,MT2,SF,true,LF,UF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct IsLower< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public BoolConstant< LF ||
                          ( IsLower_v<MT1> && IsLower_v<MT2> ) ||
                          ( ( SF || HF ) && IsUpper_v<MT1> && IsUpper_v<MT2> ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNILOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct IsUniLower< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public BoolConstant< ( IsUniLower_v<MT1> && IsUniLower_v<MT2> ) ||
                          ( ( SF || HF ) && IsUniUpper_v<MT1> && IsUniUpper_v<MT2> ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct IsStrictlyLower< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public BoolConstant< ( IsStrictlyLower_v<MT1> && IsLower_v<MT2> ) ||
                          ( IsStrictlyLower_v<MT2> && IsLower_v<MT1> ) ||
                          ( ( SF || HF ) &&
                            ( ( IsStrictlyUpper_v<MT1> && IsUpper_v<MT2> ) ||
                              ( IsStrictlyUpper_v<MT2> && IsUpper_v<MT1> ) ) ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct IsUpper< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public BoolConstant< UF ||
                          ( IsUpper_v<MT1> && IsUpper_v<MT2> ) ||
                          ( ( SF || HF ) && IsLower_v<MT1> && IsLower_v<MT2> ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct IsUniUpper< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public BoolConstant< ( IsUniUpper_v<MT1> && IsUniUpper_v<MT2> ) ||
                          ( ( SF || HF ) && IsUniLower_v<MT1> && IsUniLower_v<MT2> ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct IsStrictlyUpper< TDMatSMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public BoolConstant< ( IsStrictlyUpper_v<MT1> && IsUpper_v<MT2> ) ||
                          ( IsStrictlyUpper_v<MT2> && IsUpper_v<MT1> ) ||
                          ( ( SF || HF ) &&
                            ( ( IsStrictlyLower_v<MT1> && IsLower_v<MT2> ) ||
                              ( IsStrictlyLower_v<MT2> && IsLower_v<MT1> ) ) ) >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
