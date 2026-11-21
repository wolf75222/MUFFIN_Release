//
//  ModelWarnings.hpp
//  Muffin
//
//  Warning lines of the physical models. Outside hold() and release() they are printed by py::print as in the
//  original code; in between they are kept and release() writes them to sys.stdout with one call, the same text.
//

#ifndef ModelWarnings_hpp
#define ModelWarnings_hpp

#include <string>
#include <pybind11/pybind11.h>

namespace py = pybind11;

class ModelWarnings
{
public:
    static void hold() {s_hold = true;}

    static void release()
    {
        s_hold = false;
        if (s_text.empty()) return;
        std::string text;
        text.swap(s_text);
        py::module_::import("sys").attr("stdout").attr("write")(text);
    }

    // State out of the realizability domain of the 5M and 6M HyQMOM models. py::print joins its arguments with a
    // space and writes str() of a double, float.__repr__, which is PyOS_double_to_string(x, 'r', 0, ADD_DOT_0).
    static void realizability(double qStar, double rStar)
    {
        if (!s_hold) {
            py::print("WARNING: You are out of the realizability");
            py::print("QStar = ", qStar, "  RStar = ", rStar);
            return;
        }
        s_text += "WARNING: You are out of the realizability\nQStar =  ";
        appendRepr(qStar);
        s_text += "   RStar =  ";
        appendRepr(rStar);
        s_text += '\n';
    }

private:
    static void appendRepr(double x)
    {
        char* s = PyOS_double_to_string(x, 'r', 0, Py_DTSF_ADD_DOT_0, nullptr);
        if (!s) throw py::error_already_set();
        s_text += s;
        PyMem_Free(s);
    }

    static inline bool s_hold = false;
    static inline std::string s_text;
};

#endif /* ModelWarnings_hpp */
