#pragma once

#include <cassert>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <memory>
#include <limits>
#include <span>
#include <type_traits>


#define GCC_COMPILER (defined(__GNUC__) && !defined(__clang__))
#define CLANG_COMPILER (!defined(__GNUC__) && defined(__clang__))

#ifdef _OPENMP
#define OMP_PRAGMA(x) _Pragma(x)
#else
#define OMP_PRAGMA(x)
#endif

namespace Array {
    constexpr int ALIGN_BYTES = 512; // Alignment for memory allocation
    constexpr size_t all = std::numeric_limits<size_t>::max();
    constexpr size_t newaxis = std::numeric_limits<size_t>::max() - 1;

    namespace io {
        inline std::ostream& _newl(std::ostream& s) {s << '\n'; return s;}

        template <typename T>
        void ArrayExit(const T x)
        {
            std::cerr << _newl << "ERROR: " << x << "." << std::endl;
            exit(1);
        }

    }

    namespace alloc {
        enum alloc_state {unallocated=0, allocated=2, temporary=4};
        template<class T>
        std::shared_ptr<T[]> aligned_sptr(const size_t align, const size_t size)
        {
            void* ptr = nullptr;
            if (posix_memalign(&ptr, align, size) != 0) {
                io::ArrayExit("Memory allocation failed");
            }
            auto buffer = std::shared_ptr<T[]>(static_cast<T*>(ptr), free);
            return buffer;
        }
    }

    namespace indexing {
        // Calculate the strides for a given shape with no
        template<int D>
        inline std::array<size_t, D> calculate_strides(const std::array<size_t, D>& shape) {
            std::array<size_t, D> strides;
            const size_t ndim = shape.size();
            for (size_t i = 0; i < ndim; i++) {
                strides[i] = 1;
                for (size_t j = i+1; j < ndim; j++) {
                    strides[i] *= shape[j];
                }
            }
            return strides;
        }

        template<int D>
        inline size_t get_size_from_shape(const std::array<size_t,D>& shape) {
            size_t size = 1;
            for (int i = 0; i < D; i++) {
                size *= shape[i];
            }
            return size;
        }

    }



    template<class T, int D>
    class Array {
    private:
        const size_t size;
        const std::array<size_t,D> shape, strides;
        std::shared_ptr<T[]> buffer;
        mutable alloc::alloc_state state = alloc::alloc_state::unallocated;
        bool is_view = false;

    public:
        Array(): size(0), shape({}), buffer(nullptr) {}; // Default constructor

        // Constructor for a new array
        template<class... Dims>
        explicit Array(Dims... dims):
            size((dims * ...)),
            shape({static_cast<size_t>(dims)...}),
            strides(indexing::calculate_strides<D>(shape)),
            buffer(nullptr)
        {
            constexpr size_t ndim = sizeof...(dims);
            static_assert(ndim == D, "Incorrect size for initialization of array!");
            for (size_t i = 0; i < D; i++) {
                if (shape[i] == 0) {
                    io::ArrayExit("Array dimensions must be positive, got " + std::to_string(shape[i]) + " for dimension " + std::to_string(i) + ".");
                }
                if (shape[i] >= newaxis) {
                    io::ArrayExit("Array dimension "+ std::to_string(shape[i])+ " must be smaller than " + std::to_string(newaxis) + " for dimension" + std::to_string(i) + ".");
                }
            }

            Allocate();
        }
        // Constructor from a C-style array
        Array(T* A, std::array<size_t,D> shape) :
            size(indexing::get_size_from_shape<D>(shape)),
            shape(shape),
            strides(indexing::calculate_strides<D>(shape)),
            buffer(nullptr) {

            for (size_t i = 0; i < D; i++) {
                if (shape[i] == 0) {
                    io::ArrayExit("Array dimensions must be positive, got " + std::to_string(shape[i]) + " for dimension " + std::to_string(i) + ".");
                }
                if (shape[i] >= newaxis) {
                    io::ArrayExit("Array dimension "+ std::to_string(shape[i])+ " must be smaller than " + std::to_string(newaxis) + " for dimension" + std::to_string(i) + ".");
                }
            }
            Allocate();
            this->Load(A);
        }
        // Constructor from an Array - makes a copy
        Array(const Array<T,D>& A) :
            size(A.Size()),
            shape(A.Shape()),
            strides(indexing::calculate_strides<D>(A.Shape())),
            buffer(nullptr) {

            Allocate();
            this->Load(A.Buffer());
        }

        // Copy and move assignment
        Array& operator = (T *a) {this->Load(a); return *this;}
        Array& operator = (const Array& A) {
            CheckShape(this,A);
            this->Load(A.Buffer());
            return *this;
        }
        Array& operator = (Array&& A) noexcept {
            CheckShape(this,A);
            this->buffer = std::move(A.Buffer());
            A.Deallocate();
            return *this;
        }

        template <size_t array_size>
        Array& View (std::span<size_t, array_size> dims) {
            if (array_size >= D) {
                io::ArrayExit("View dimensions must be less than the number of dimensions of the array");
            }
            Array<T, D - array_size> output;


            return *output;
        }

        Array&  View (std::initializer_list<size_t> dims) {
            View(std::span{dims.begin(), dims.size()});
        }

        template<typename BufferType>
        void Load(BufferType buffer) {
            OMP_PRAGMA("omp parallel for default(none) shared(buffer)")
            for(size_t i = 0; i < this->size; i++) {
                this->buffer[i] = buffer[i];
            }
        }


        [[nodiscard]] const alloc::alloc_state& State() const {
            return state;
        }
        [[nodiscard]] const bool& IsView() const {
            return is_view;
        }
        [[nodiscard]] std::shared_ptr<T[]> Buffer() const {
            return buffer;
        }

        // Allocator
        void Allocate() {
            if (is_view) {
                io::ArrayExit("Cannot allocate a view");
            }
            if (state == alloc::unallocated) {
                buffer=alloc::aligned_sptr<T>(ALIGN_BYTES, size*sizeof(T));
                state = alloc::alloc_state::allocated;
            }
            else {
                io::ArrayExit("Array already allocated");
            }
        };

        void Deallocate() {
            if (is_view) {
                io::ArrayExit("Cannot allocate a view");
            }
            if (state == alloc::allocated) {
                buffer = nullptr;
                state = alloc::unallocated;
            }
            else {
                io::ArrayExit("Array already allocated");
            }
        };
        // Dimension info
        [[nodiscard]] size_t Size() const {
            return size;
        };

        [[nodiscard]] size_t Shape(const int i) const {
            return shape[i];
        }

        [[nodiscard]] const std::array<size_t,D>& Shape() const {
            return shape;
        }

        [[nodiscard]] static size_t NDim() {
            return D;
        }
        // Indexing operations
        static void CheckIndex(size_t const & index, size_t const & dimension) {
#ifdef DEBUG
            if (((index > 0 && index >= shape[dimension])) || ((index < 0) && index < shape[dimension])) {
                std::stringstream ss;
                ss << "Index " << index << " out of bounds for dimension of size " << shape[dimension];
                io::ArrayExit(ss.str().c_str());
            }
#endif
        }

        void CheckShape(const Array& other) {
#ifdef DEBUG
            bool same = other.NDim() == NDim();
            if(same)
                for (int i = 0; i < D; i++) {
                    same &= (this->Shape(i) == other.Shape(i));
                }
            if (!same) {
                io::ArrayExit("Arrays do not have the same shape!");
            }
#endif
        }


        template <typename ... DIMS>
        size_t _get_linear_index(DIMS const & ... dims) const {
            constexpr size_t ndim = sizeof...(dims) ;
            static_assert(ndim == D, "Incorrect number of indices for slicing array!");
            const std::array<long, sizeof...(DIMS)> index{dims...};
            size_t linear_index = 0;
            for (size_t i = 0; i < D; i++) {
                const size_t idx = (index[i] >= 0) ? index[i] : shape[i] + index[i];
                CheckIndex(idx, i);
                linear_index += idx * strides[i];
            }
            return linear_index;
        }

        // Return the scalar value at the given index (i,j,k,...)
        template<class... Dims>
        decltype(auto) operator[](this auto& self, Dims...dims) {
            constexpr size_t ndim = sizeof...(dims) ;
            static_assert(ndim == D, "Incorrect number of indices for slicing array!");
            const size_t index = self._get_linear_index(dims...);
            return self.buffer[index];
        }



        // Arithmetic operations

        // Addition
        template<typename Floating>
        Array& operator+= (const Floating a) {
            OMP_PRAGMA("omp parallel for default(none) shared(a)")
            for (size_t i = 0; i < size; i++) {
                buffer[i] += a;
            }
            return *this;
        }
        Array& operator+= (const Array& a) {
            CheckShape(a);
            auto a_buffer = a.Buffer().get();
            OMP_PRAGMA("omp parallel for default(none) shared(a_buffer)")
            for (size_t i = 0; i < size; i++) {
                buffer[i] += a_buffer[i];
            }
            return *this;
        }
        template<typename RHS>
        friend Array operator+ (Array lhs, const RHS& rhs) {
            lhs += rhs;
            return lhs;
        }

        // Subtraction
        template<typename Floating>
        Array& operator-= (const Floating a) {
            OMP_PRAGMA("omp parallel for default(none) shared(a)")
            for (size_t i = 0; i < size; i++) {
                buffer[i] -= a;
            }
            return *this;
        }
        Array& operator-= (const Array& a) {
            CheckShape(a);
            auto a_buffer = a.Buffer().get();
            OMP_PRAGMA("omp parallel for default(none) shared(a_buffer)")
            for (size_t i = 0; i < size; i++) {
                buffer[i] -= a_buffer[i];
            }
            return *this;
        }
        template<typename RHS>
        friend Array operator- (Array lhs, const RHS& rhs) {
            lhs -= rhs;
            return lhs;
        }

        // Mult
        template<typename Floating>
        Array& operator*= (const Floating a) {
            OMP_PRAGMA("omp parallel for default(none) shared(a)")
            for (size_t i = 0; i < size; i++) {
                buffer[i] *= a;
            }
            return *this;
        }
        Array& operator*= (const Array& a) {
            CheckShape(a);
            auto a_buffer = a.Buffer().get();
            OMP_PRAGMA("omp parallel for default(none) shared(a_buffer)")
            for (size_t i = 0; i < size; i++) {
                buffer[i] *= a_buffer[i];
            }
            return *this;
        }
        template<typename RHS>
        friend Array operator* (Array lhs, const RHS& rhs) {
            lhs *= rhs;
            return lhs;
        }

        // Div
        template<typename Floating>
        Array& operator/= (const Floating a) {
            OMP_PRAGMA("omp parallel for default(none) shared(a)")
            for (size_t i = 0; i < size; i++) {
                buffer[i] /= a;
            }
            return *this;
        }
        Array& operator/= (const Array& a) {
            CheckShape(a);
            auto a_buffer = a.Buffer().get();
            OMP_PRAGMA("omp parallel for default(none) shared(a_buffer)")
            for (size_t i = 0; i < size; i++) {
                buffer[i] /= a_buffer[i];
            }
            return *this;
        }
        template<typename RHS>
        friend Array operator/ (Array lhs, const RHS& rhs) {
            lhs /= rhs;
            return lhs;
        }
        // Reductions
        double Sum() {
            double reduction = 0.0;
            OMP_PRAGMA("omp parallel for default(none) reduction(+:reduction)")
            for (size_t i = 0; i < size; i++) {
                reduction += buffer[i];
            }
            return reduction;
        };
        double Prod() {
            double reduction = 1.0;
            OMP_PRAGMA("omp parallel for default(none) reduction(*:reduction)")
            for (size_t i = 0; i < size; i++) {
                reduction *= buffer[i];
            }
            return reduction;
        };
        T Min() {
            T reduction = std::numeric_limits<T>::max();
            OMP_PRAGMA("omp parallel for default(none) reduction(min:reduction)")
            for (size_t i = 0; i < size; i++) {
                if (buffer[i] < reduction) reduction = buffer[i];
            }
            return reduction;
        };
        T Max() {
            T reduction = std::numeric_limits<T>::min();
            OMP_PRAGMA("omp parallel for default(none) reduction(max:reduction)")
            for (size_t i = 0; i < size; i++) {
                if (buffer[i] > reduction) reduction = buffer[i];
            }
            return reduction;
        };
        double Moment(const int m) {
            double reduction = 0.0;
            OMP_PRAGMA("omp parallel for default(none) reduction(+:reduction)")
            for (size_t i = 0; i < size; i++) {
                reduction += std::pow(buffer[i], m);
            }
            return reduction / size;
        };
        double Mean() {
            return Moment(1);
        };
        double Variance() {
            return Moment(2) - std::pow(Moment(1), 2);
        }
        double Std() {
            return std::sqrt(Variance());
        }

    };


}