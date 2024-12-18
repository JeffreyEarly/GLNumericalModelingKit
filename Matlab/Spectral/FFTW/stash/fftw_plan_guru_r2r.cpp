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

    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    std::ostringstream stream;
    
    void flattenDimensionsAndTransformIndices(std::vector<size_t>& indims, std::vector<size_t>& transformDims) {
        for (auto& num : transformDims) {
            num -= 1; // Convert 1-based to 0-based indexing
        }

        std::vector<size_t> dims;
        for (size_t i = 0; i < indims.size(); ++i) {
            if (indims[i] != 1) {
                dims.push_back(indims[i]);
            } else {
                for (auto& tDim : transformDims) {
                    if (tDim == i) throw std::invalid_argument("Transform dimension has length 1.");
                    if (tDim > i) --tDim;
                }
            }
        }
        indims = dims;
    }

    void flattenDimensions(std::vector<size_t>& indims) {
        indims.erase(std::remove(indims.begin(), indims.end(), 1), indims.end());
    }
    
    void fftwDimsSetup(const std::vector<size_t>& dims, const std::vector<size_t>& transformDims,
                       std::vector<fftw_iodim>& iodims, std::vector<fftw_iodim>& howmany_dims) {
        std::vector<size_t> outdims(dims);

        std::vector<fftw_iodim> alldims;
        for (size_t i = 0; i < dims.size(); i++) {
            fftw_iodim newdim;
            newdim.n = dims[i];
            if (alldims.size() == 0) {
                newdim.is = 1;
                newdim.os = 1;
            } else {
                size_t n = alldims.size();
                newdim.is = alldims[n-1].is*dims[n-1];
                newdim.os = alldims[n-1].os*outdims[n-1];
            }
            alldims.push_back(newdim);
        }

        iodims.clear();
        howmany_dims.clear();

        for (size_t i = 0, t = 0; i < alldims.size(); ++i) {
            if (t < transformDims.size() && i == transformDims[t]) {
                iodims.push_back(alldims[i]);
                ++t;
            } else {
                howmany_dims.push_back(alldims[i]);
            }
        }
    }
    
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
        
        auto realMatrixDims = typedArrayToSizeTVector(inputs[0]);
        auto transformDims = typedArrayToSizeTVector(inputs[1]);
        std::vector<size_t> flatRealMatrixDims = realMatrixDims;
        flattenDimensionsAndTransformIndices(flatRealMatrixDims, transformDims);
        
        std::vector<fftw_iodim> iodims, howmany_dims;
        fftwDimsSetup(flatRealMatrixDims, transformDims, iodims, howmany_dims);
        
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
        
        printDims(iodims);
        printDims(howmany_dims);
        
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

    void printDims(const std::vector<fftw_iodim>& alldims) {
        matlab::data::ArrayFactory factory;
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        for (size_t i = 0; i < alldims.size(); i++) {
            matlabPtr->feval(u"fprintf", 0,
                     std::vector<matlab::data::Array>(
                      {factory.createScalar("alldims (%d, %d, %d)\n"),
                          factory.createScalar(alldims[i].n),
                          factory.createScalar(alldims[i].is),
                          factory.createScalar(alldims[i].os)}));
        }
    }

};



