//
//  create_dct_plan_mex.cpp
//  
//
//  Created by Jeffrey Early on 12/7/24.
//

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "fftw3.h"
//#include "omp.h"
#include <stdexcept>
#include <vector>

using namespace matlab::mex;
using namespace matlab::data;

/* Helper function to generate an error message from given string,
 * and display it over MATLAB command prompt.
 */
class MexFunction : public matlab::mex::Function {
private:
    struct PlanHandle {
        fftw_plan plan;
    };

    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr;
    std::ostringstream stream;
    
public:
    /* Constructor for the class. */
    MexFunction()
    {
      matlabPtr = getEngine();
    }
    
    /* Helper function to generate an error message from given string,
     * and display it over MATLAB command prompt.
     */
    void displayError(std::string errorMessage)
    {
        ArrayFactory factory;
        matlabPtr->feval(u"error", 0, std::vector<Array>({factory.createScalar(errorMessage) }));
    }
    
    /* Helper function to print output string on MATLAB command prompt. */
//    void displayOnMATLAB(std::ostringstream stream)
//    {
//      ArrayFactory factory;
//      matlabPtr->feval(u"fprintf", 0, std::vector<Array>({ factory.createScalar(stream.str())}));
//    }
    void displayOnMATLAB(std::ostringstream& stream) {
            // Pass stream content to MATLAB fprintf function
        ArrayFactory factory;
            matlabPtr->feval(u"fprintf", 0, std::vector<Array>({ factory.createScalar(stream.str()) }));
            // Clear stream buffer
            stream.str("");
        }
    
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) override {
        namespace mat = matlab::data;
        mat::ArrayFactory factory;
        
        // Validate input arguments
        if (inputs.size() != 5 || outputs.size() != 1) {
            displayError("Usage: plan = create_dct_plan_mex(rank, dims, howmany_rank)");
            return;
        }
        stream << "Here 1\n";
        displayOnMATLAB(stream);
        // Extract inputs
        int nCores = static_cast<int>(inputs[2][0]);
        fftw_r2r_kind kind = (fftw_r2r_kind) static_cast<int>(inputs[3][0]);
        unsigned planner = static_cast<unsigned>(inputs[4][0]);
  
        stream << "Here 2\n";
        displayOnMATLAB(stream);
        
        // Convert dimsStructArray to FFTW's fftw_iodim structure
        const mat::StructArray dimsStructArray = inputs[0];
        int rank = static_cast<int>(dimsStructArray.getNumberOfElements());
        std::vector<fftw_iodim> dims(rank);
        for (int i = 0; i < rank; ++i) {
            matlab::data::TypedArray<double> field1 = dimsStructArray[i]["n"];
            matlab::data::TypedArray<double> field2 = dimsStructArray[i]["is"];
            matlab::data::TypedArray<double> field3 = dimsStructArray[i]["os"];
            dims[i].n = (int) field1[0];
            dims[i].is = (int) field2[0];
            dims[i].os = (int) field3[0];
        }
        stream << "Here 3\n";
        displayOnMATLAB(stream);
        
        int howmany_rank = 0;
        std::vector<fftw_iodim> howmany_dims(0);
        if (inputs[1].getNumberOfElements() != 0) {
            stream << "Here 3.5\n";
            displayOnMATLAB(stream);
            const mat::StructArray howmany_dimsStructArray = inputs[1];
            howmany_rank = static_cast<int>(howmany_dimsStructArray.getNumberOfElements());
            howmany_dims.resize(howmany_rank);
            for (int i = 0; i < rank; ++i) {
                matlab::data::TypedArray<double> field1 = howmany_dimsStructArray[i]["n"];
                matlab::data::TypedArray<double> field2 = howmany_dimsStructArray[i]["is"];
                matlab::data::TypedArray<double> field3 = howmany_dimsStructArray[i]["os"];
                howmany_dims[i].n = (int) field1[0];
                howmany_dims[i].is = (int) field2[0];
                howmany_dims[i].os = (int) field3[0];
            }
        }
        
        stream << "Here 4\n";
        displayOnMATLAB(stream);
        // Allocate input and output arrays
        int totalSize = 1;
        for (const auto& dim : dims) {
            totalSize *= dim.n;
            stream << "dim size " << dim.n << "\n";
            displayOnMATLAB(stream);
        }
        for (const auto& dim : howmany_dims) {
            totalSize *= dim.n;
            stream << "dim size " << dim.n << "\n";
            displayOnMATLAB(stream);
        }
        stream << "Here 5 and allocated size of " << totalSize << "\n";
        displayOnMATLAB(stream);
        
//        omp_set_num_threads(nCores);
        fftw_init_threads();
        fftw_plan_with_nthreads(nCores);
        
        double* in = fftw_alloc_real(totalSize);
        double* out = fftw_alloc_real(totalSize);
        fftw_plan plan = fftw_plan_guru_r2r(rank, dims.data(), howmany_rank, howmany_dims.data(), in, out, &kind, planner);
        fftw_free(in);
        fftw_free(out);
        
        stream << "Here 6\n";
        displayOnMATLAB(stream);
        if (!plan) {
            displayError("Failed to create FFTW plan.");
            return;
        }

        // Wrap the plan and buffers in a struct
        PlanHandle* handle = new PlanHandle{plan};

        // Return the handle as an opaque pointer
        outputs[0] = factory.createScalar(reinterpret_cast<uint64_t>(handle));
    }


};



