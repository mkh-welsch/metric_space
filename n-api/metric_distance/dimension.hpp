#ifndef DIMENSION_HPP
#define DIMENSION_HPP

#include <vector>
#include <variant> // for DimensionSet

#include "metric_distance.hpp"



namespace metric
{



    template <
        typename Record, // needed for accessors
        //typename Container, // type of field, that metric functor takes as input
        //template <typename> class Metric // metric functor template
        typename MetricType
        >
    class Dimension
    {
    public:
        using InputValueType = typename MetricType::value_type;
        using ReturnValueType = typename MetricType::distance_type;

        MetricType DistanceFunctor; // = MetricType();

        template<typename RetType, typename RecordType, typename M>
        Dimension(std::vector<std::function<RetType(RecordType)>> accessors_,  const M & m_):DistanceFunctor(m_) // vector of accessors as input
        {
            accessors = accessors_;
        }

        std::vector<InputValueType> operator ()(Record record)
        {
            std::vector<InputValueType> result;
            for (size_t i=0; i<accessors.size(); i++)
            {
                result.push_back(accessors[i](record));
            }
            return result;
        }

        ReturnValueType get_distance(Record r1, Record r2)
        {
            std::vector<InputValueType> field1 = (*this)(r1);
            std::vector<InputValueType> field2 = (*this)(r2);
            // auto d_dist_vis = [&, this, field1, field2]()
            //                   {
            //                       return DistanceFunctor(field1, field2); // here we support only 2 operands
            //                   };
            // // TODO add static_assert
            // ReturnValueType distance = std::visit(d_dist_vis);
//            return distance;
            return DistanceFunctor(field1, field2);
        }

    private:
        std::vector<std::function<InputValueType(Record)>> accessors; // first by field, than by elements of field
    };
    template<typename R, typename Rt, typename M>
    Dimension(std::vector<std::function<R(Rt)>> accessors_,  const M & m_) -> Dimension<Rt,M>;
// // code for DimensionSet

//template <typename Record, typename Container, typename DistType>
//class DimensionSet
//{

//private:

//    typedef typename Container::value_type InputValueType; // elementary values inside field we create

//    typedef  std::variant<
//      metric::Dimension<Record, std::vector<InputValueType>, metric::distance::Euclidian>,
//      metric::Dimension<Record, std::vector<InputValueType>, metric::distance::Manhatten>,
//      metric::Dimension<Record, std::vector<InputValueType>, metric::distance::P_norm>,
//      metric::Dimension<Record, std::vector<InputValueType>, metric::distance::Euclidian_thresholded>,
//      metric::Dimension<Record, std::vector<InputValueType>, metric::distance::Cosine>,
//      //metric::Dimension<Record, std::vector<InputValueType>, metric::distance::SSIM>,
//      //metric::Dimension<Record, std::vector<InputValueType>, metric::distance::TWED>,
//      //metric::Dimension<Record, std::vector<InputValueType>, metric::distance::EMD>, // TODO determine why return type deduction fails
//      metric::Dimension<Record, std::vector<InputValueType>, metric::distance::Edit>
//    > VariantType;

//    std::vector<VariantType> variant_dims;

//public:

//    typedef Container FieldType;
//    typedef DistType DistanceType;

//    template <template <typename> class Metric>
//    void add(Dimension<Record, Container, Metric> dim) // add Dimension object to the set
//    {
//        variant_dims.push_back(dim);
//    }

//    size_t size()
//    {
//        return variant_dims.size();
//    }

//    Container get(Record r, size_t dim_idx) // get dim_idx field
//    {
//        Container field;
//        if (dim_idx >= variant_dims.size() || dim_idx < 0)
//            return field; // of rise exception?
//        auto d_vis = [&r](auto & d) { return d(r); };
//        field = std::visit(d_vis, variant_dims[dim_idx]); // call Dimension functor for current Record
//        return field;
//    }

//    DistanceType get_distance(Record r1, Record r2, size_t dim_idx)
//    {
//        Container field1 = get(r1, dim_idx);
//        Container field2 = get(r2, dim_idx);
//        auto d_dist_vis = [field1, field2](auto & d)
//        {
//            return (DistanceType)d.DistanceFunctor(field1, field2);
//        };
//        // TODO add static_assert
//        DistanceType distance = std::visit(d_dist_vis, variant_dims[dim_idx]);
//        return distance;
//    }
//};





} // namespace metric


#endif // DIMENSION_HPP
