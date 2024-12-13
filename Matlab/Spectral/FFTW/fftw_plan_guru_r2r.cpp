//
//  create_dct_plan_mex.cpp
//  
//
//  Created by Jeffrey Early on 12/7/24.
//

#include "mex.hpp"
#include "mexAdapter.hpp"
#include "fftw3.h"
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
    
    // Function to compute FFTW iodim structure based on a vector of dimensions
    void createFFTWIODim(const std::vector<size_t>& dims, int transformDim, std::vector<fftw_iodim>& iodims, std::vector<fftw_iodim>& howmany_dims) {
        // Validate transformDim
        if (transformDim < 1 || transformDim > static_cast<int>(dims.size())) {
            throw std::invalid_argument("Invalid transform dimension");
        }
        
        // Convert 1-based indexing to 0-based indexing
        transformDim -= 1;
        
        // This loop computes all the strides for the matrix and also eliminates any singleton dimensions.
        std::vector<fftw_iodim> alldims;
        for (size_t i = 0; i < dims.size(); i++) {
            if (dims[i] == 1 && transformDim >= i) {
                transformDim -= 1;
            } else if (dims[i] == 1 && transformDim == i) {
                throw std::invalid_argument("Transform dimension has length 1.");
            }
            if (dims[i] == 1) {
                continue;
            }
            fftw_iodim newdim;
            newdim.n = static_cast<int>(dims[i]);
            if (alldims.size() == 0) {
                newdim.is = 1;
                newdim.os = 1;
            } else {
                size_t n = alldims.size();
                newdim.is = alldims[n].is*alldims[n].n;
                newdim.os = alldims[n].os*alldims[n].n;
            }
            alldims.push_back(newdim);
        }

        iodims.clear();
        howmany_dims.clear();
        for (size_t i = 0; i < alldims.size(); ++i) {
            if (i == transformDim) {
                iodims.push_back(alldims[i]);
            } else {
                howmany_dims.push_back(alldims[i]);
            }
        }
    }

    std::vector<size_t> typedArrayToSizeTVector(const matlab::data::TypedArray<double>& inputArray) {
        std::vector<size_t> dims(inputArray.getNumberOfElements());
        size_t index = 0;
        for (const auto& val : inputArray) {
            dims[index++] = static_cast<size_t>(val);
        }
        return dims;
    }
    
    std::vector<size_t> arrayToSizeTVector(const matlab::data::Array& inputArray) {
        if (inputArray.getType() != matlab::data::ArrayType::DOUBLE) {
            throw std::invalid_argument("Input array must be of type double.");
        }

        const matlab::data::TypedArray<double>& typedArray = static_cast<const matlab::data::TypedArray<double>&>(inputArray);
        std::vector<size_t> dims(typedArray.getNumberOfElements());

        size_t index = 0;
        for (const auto& val : typedArray) {
            if (val < 0) {
                throw std::invalid_argument("Dimensions must be non-negative.");
            }
            dims[index++] = static_cast<size_t>(val);
        }

        return dims;
    }

    
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
        
        std::vector<size_t> dims = typedArrayToSizeTVector(inputs[0]);
        size_t transformDim = static_cast<size_t>(inputs[1][0]);
        std::vector<fftw_iodim> iodims;
        std::vector<fftw_iodim> howmany_dims;
        createFFTWIODim(dims, transformDim, iodims, howmany_dims);
        
        // Extract inputs
        int nCores = static_cast<int>(inputs[2][0]);
        fftw_r2r_kind kind = (fftw_r2r_kind) static_cast<int>(inputs[3][0]);
        unsigned planner = static_cast<unsigned>(inputs[4][0]);
        
        stream << "Here 4\n";
        displayOnMATLAB(stream);
        // Allocate input and output arrays
        int totalSize = 1;
        for (const auto& dim : iodims) {
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
        
        fftw_init_threads();
        fftw_plan_with_nthreads(nCores);
        
        double* in = fftw_alloc_real(totalSize);
        double* out = fftw_alloc_real(totalSize);
        fftw_plan plan = fftw_plan_guru_r2r(static_cast<int>(iodims.size()), iodims.data(), static_cast<int>(howmany_dims.size()), howmany_dims.data(), in, out, &kind, planner);
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



