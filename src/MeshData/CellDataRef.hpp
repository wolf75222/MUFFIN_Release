//
//  CellDataRef.hpp
//  Muffin
//
//  Created by Alejandro Alvarez Laguna on 15/08/18.
//  Copyright © 2018 Alejandro Alvarez Laguna. All rights reserved.
//


#ifndef CellDataRef_hpp
#define CellDataRef_hpp

#include <cstddef>
#include <stdio.h>
#include "../Parameters.hpp"
#include <pybind11/numpy.h>
namespace py = pybind11;

// Column iCell of a 2D array (equations x cells). The address of the column and the stride between two
// equations are read once from the array, whose buffer is allocated once in Simulation1D::setup and never
// moved: operator[] is a plain load, the same value as m_array->unchecked<2>()(i + offset, iCell).
class CellDataRef
{
public:
    CellDataRef() : m_array(nullptr), m_iCell(0), m_offset(0), m_ptr(nullptr), m_stride(0) {}
    CellDataRef(py::array_t<double>* array, int iCell, int offset = 0) : m_array(array), m_iCell(iCell), m_offset(offset)
    {
        m_stride = array->strides(0) / static_cast<py::ssize_t>(sizeof(double));
        const py::ssize_t cell = array->strides(1) / static_cast<py::ssize_t>(sizeof(double));
        m_ptr = array->mutable_data() + offset * m_stride + iCell * cell;
    }
    CellDataRef withOffset(int offset) const {return CellDataRef(m_array, m_iCell,offset);}
    py::array_t<double> getView() const {
        py::buffer_info info = m_array->request();
        size_t itemSize = info.itemsize;
        size_t numRows = info.shape[0];
        size_t numCols = info.shape[1];
        size_t stride = info.strides[1];

        double* dataPtr = static_cast<double*>(info.ptr);
        double* columnPtr = dataPtr + m_iCell * stride / itemSize;

        return py::array_t<double>(
            {numRows},
            {numCols * itemSize},
            columnPtr
        );
    }
    double operator [](int i) const {return m_ptr[i*m_stride];}
    double & operator [](int i) {return m_ptr[i*m_stride];}
    const double* data() const {return m_ptr;}
    std::ptrdiff_t stride() const {return m_stride;}

protected:
    py::array_t<double>* m_array;
    int m_iCell;
    int m_offset;
    double* m_ptr;
    std::ptrdiff_t m_stride;
};

#endif /* CellDataRef_hpp */
