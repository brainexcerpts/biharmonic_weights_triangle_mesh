#include "solvers.hpp"

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include <Eigen/QR>
#include <Eigen/LU>
#include <Eigen/SparseLU>

// -----------------------------------------------------------------------------

typedef Eigen::Triplet<double, int> Triplet;
/// declares a column-major sparse matrix type of double
typedef Eigen::SparseMatrix<double> Sparse_mat;

// -----------------------------------------------------------------------------

/// Are angles obtuse between e0 and e1 (meaning a < PI/2)
static bool check_obtuse(const Vec3& e0,
                         const Vec3& e1)
{
    return e0.normalized().dot( e1.normalized() ) >= 0.f;
}

// -----------------------------------------------------------------------------

/**
 * @brief Check if a triangle is obtuse (every angles < PI/2)
 * Here how edges must be defined :
   @code
      p0
      |\
      | \
      |  \
      |   \
      |____\
     p1     p2

    Vec3 e0 = p1 - p0;
    Vec3 e1 = p2 - p0;
    Vec3 e2 = p2 - p1;
   @endcode
 */
static bool check_obtuse(const Vec3& e0,
                         const Vec3& e1,
                         const Vec3& e2)
{
    return check_obtuse(e0, e1) && check_obtuse(-e0, e2) && check_obtuse(-e1, -e2);
}

// -----------------------------------------------------------------------------

inline static
double mixed_voronoi_area(const Vec3& pi,
                          const Vec3& pj0,
                          const Vec3& pj1)
{
    double area = 0.;
    Vec3 e0 = pj0 - pi ;
    Vec3 e1 = pj1 - pi;
    Vec3 e2 = pj1 - pj0;

    if( check_obtuse(e0, e1, e2) )
    {
        area = (1/8.) * (double)(e0.norm_squared() * (-e0).cotan( e2) +
                                 e1.norm_squared() * (-e1).cotan(-e2));
    }
    else
    {
        const double ta = (double)e0.cross( e1 ).norm() / 2.; // Tri area
        area = ta / (check_obtuse(e0, e1) ? 2. : 4.);
    }
    return area;
}

// -----------------------------------------------------------------------------

static
double get_cell_area(int vidx,
                     const std::vector< Vec3 >& vertices,
                     const std::vector< std::vector<int> >& edges )
{
    double area = 0.0;
    const Vec3 c_pos = vertices[vidx];
    //get triangles areas
    for(int e = 0; e < (int)edges[vidx].size(); ++e)
    {
        int ne = (e + 1) % edges[vidx].size();

        if(false){
            Vec3 p0 = vertices[edges[vidx][e ] ];
            Vec3 p1 = vertices[edges[vidx][ne]];
            area += mixed_voronoi_area(c_pos, p0, p1 );
        }else{
            Vec3 edge0 = vertices[edges[vidx][e ]] - c_pos;
            Vec3 edge1 = vertices[edges[vidx][ne]] - c_pos;
            area += (edge0.cross(edge1)).norm() / 3.f;
        }
    }
    return area;
}


// -----------------------------------------------------------------------------

/// @return A sparse representation of the Laplacian
/// list[ith_row][list of columns] = Triplet(ith_row, jth_column, matrix value)
static
std::vector<std::vector<Triplet>>
get_laplacian(const std::vector< Vec3 >& vertices,
              const std::vector< std::vector<int> >& edges )
{
    std::cout << "BUILD LAPLACIAN MATRIX" << std::endl;
    unsigned nv = unsigned(vertices.size());
    std::vector<std::vector<Triplet>> mat_elemts(nv);
    for(int i = 0; i < nv; ++i)
        mat_elemts[i].reserve(10);

    for(int i = 0; i < nv; ++i)
    {
        const Vec3 c_pos = vertices[i];
        //double area = get_cell_area(i, vertices, edges);

        //area = 1. / ((1e-10 + area) * 2.f);

        //get laplacian
        double sum = 0.;
        int nb_edges = edges[i].size();
        for(int e = 0; e < nb_edges; ++e)
        {
            int next_edge = (e + 1           ) % nb_edges;
            int prev_edge = (e + nb_edges - 1) % nb_edges;

            Vec3 v1 = c_pos                 - vertices[edges[i][prev_edge]];
            Vec3 v2 = vertices[edges[i][e]] - vertices[edges[i][prev_edge]];
            Vec3 v3 = c_pos                 - vertices[edges[i][next_edge]];
            Vec3 v4 = vertices[edges[i][e]] - vertices[edges[i][next_edge]];

            double cotan1 = (v1.dot(v2)) / (1e-6 + (v1.cross(v2)).norm() );
            double cotan2 = (v3.dot(v4)) / (1e-6 + (v3.cross(v4)).norm() );

            double w = (cotan1 + cotan2)*0.5;
            //w = 1.0f;/////////////////////////////////////////////////////

            sum += w;

            mat_elemts[i].push_back( Triplet(i, edges[i][e], w) );
        }

        mat_elemts[i].push_back( Triplet(i, i, -sum) );
    }
    return mat_elemts;
}

//------------------------------------------------------------------------------

std::vector<Triplet> get_mass_matrix(const std::vector< Vec3 >& vertices,
                                     const std::vector< std::vector<int> >& edges)
{
    std::cout << "BUILD MASS MATRIX" << std::endl;
    unsigned nv = unsigned(vertices.size());
    std::vector<Triplet> triplets(nv);
    for(unsigned i = 0; i < nv; ++i) {
        double area = get_cell_area(i, vertices, edges);
        area = 1. / ((1e-10 + area));
        triplets[i] = Triplet(i, i, area);
    }
    return triplets;
}

//------------------------------------------------------------------------------

// Compute biharmonic weights
void solve_bilaplacian(const std::vector< Vec3 >& vertices,
                       const std::vector< std::vector<int> >& edges,
                       const std::vector<std::pair<Vert_idx, float> >& boundaries,
                       std::vector<double>& harmonic_weight_map)
{

    std::cout << "COMPUTE BI-LAPLACIAN WEIGHT MAP" << std::endl;

    int nv = vertices.size();
    // compute laplacian matrix of the mesh

    std::vector<std::vector<Triplet>> mat_elemts = get_laplacian(vertices, edges);



#if 0
    Eigen::MatrixXd L = Eigen::MatrixXd::Constant(nv, nv, 0.);
    for( const std::vector<Triplet>& row : mat_elemts)
        for( const Triplet& elt : row )
            L(elt.row(), elt.col()) = elt.value();

    //Eigen::ColPivHouseholderQR<Eigen::MatrixXd> llt;
    Eigen::FullPivLU<Eigen::MatrixXd> solver;
    std::cout << "BEGIN MATRIX FACTORIZATION" << std::endl;
    solver.compute( L );
    std::cout << "END MATRIX FACTORIZATION" << std::endl;

#else
    Eigen::SparseMatrix<double, Eigen::RowMajor> L(nv, nv);
    // Convert to triplets
    std::vector<Triplet> triplets;
    triplets.reserve(nv * 10);
    for( const std::vector<Triplet>& row : mat_elemts)
        for( const Triplet& elt : row )
            triplets.push_back( elt );

    L.setFromTriplets(triplets.begin(), triplets.end());

    // TODO : try L*M*L by decoupling areas into a mass matrix M.
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(nv, nv);
    triplets = get_mass_matrix(vertices, edges);
    M.setFromTriplets(triplets.begin(), triplets.end());
// Because we solve for the laplace equation we can ignore cell areas


    L = L*M*L; //////////////////////////////////////////


    // Set boundary conditions
    Eigen::VectorXd rhs = Eigen::VectorXd::Constant(nv, 0.);
    // Initialize handle
    for(const std::pair<int, float>& elt : boundaries){
        rhs( elt.first ) = double(elt.second);
        // Set row to 0.0f
        for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(L, elt.first); it; ++it) {
            it.valueRef() = 0.0;
        }
        // Set
        L.coeffRef(elt.first, elt.first) = 1.0;
    }
    Eigen::SparseMatrix<double, Eigen::ColMajor> mat = L;
    mat.makeCompressed();

    Eigen::SparseLU<Sparse_mat> solver;
    std::cout << "BEGIN SPARSE MATRIX FACTORIZATION" << std::endl;
    solver.compute( mat );
    std::cout << "END SPARSE MATRIX FACTORIZATION" << std::endl;
#endif

    harmonic_weight_map.resize(nv);
    Eigen::VectorXd res = solver.solve( rhs );
    for(int i = 0; i < nv; ++i)
        harmonic_weight_map[i] = res(i);

    return;
}

// -----------------------------------------------------------------------------

void solve_bilaplacian_diffusion(
        const std::vector< Vec3 >& vertices,
        const std::vector< std::vector<int> >& edges,
        const std::vector<std::pair<Vert_idx, float> >& boundaries,
        std::vector<double>& harmonic_weight_map)
{
    std::cout << "COMPUTE BI-LAPLACIAN DIFFUSION" << std::endl;

    int nv = vertices.size();
    // compute laplacian matrix of the mesh

    std::vector<std::vector<Triplet>> mat_elemts = get_laplacian(vertices, edges);


    Eigen::SparseMatrix<double, Eigen::RowMajor> L(nv, nv);
    // Convert to triplets
    std::vector<Triplet> triplets;
    triplets.reserve(nv * 10);
    for( const std::vector<Triplet>& row : mat_elemts)
        for( const Triplet& elt : row )
            triplets.push_back( elt );

    L.setFromTriplets(triplets.begin(), triplets.end());

    // TODO : try L*M*L by decoupling areas into a mass matrix M.
    Eigen::SparseMatrix<double, Eigen::RowMajor> M(nv, nv);
    triplets = get_mass_matrix(vertices, edges);
    M.setFromTriplets(triplets.begin(), triplets.end());

    // Because we solve for the laplace equation we can ignore cell areas

    L = L*M*L; //////////////////////////////////////////


    harmonic_weight_map.assign(nv, 0.0f);
    // Set boundary conditions
    //Eigen::VectorXd rhs = Eigen::VectorXd::Constant(nv, 0.);
    // Initialize handle
    for(const std::pair<int, float>& elt : boundaries){
        //rhs( elt.first ) = double(elt.second);
        harmonic_weight_map[elt.first] = elt.second;
        // Set row to 0.0f
        for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(L, elt.first); it; ++it) {
            it.valueRef() = 0.0;
        }
        // Set
        L.coeffRef(elt.first, elt.first) = 1.0;
    }
    //std::cout << Eigen::MatrixXd(L) << std::endl;

    std::cout << "nb verts: " << nv << std::endl;
    std::cout << "outer size" << L.outerSize() << std::endl;
    double t = 0.60;//biharmonic: 0.6f tri: 0.4
    std::vector<double> buffer = harmonic_weight_map;
    unsigned iter = 0;
    for(; iter < 10000; ++iter)
    {
        // Look up rows
        for(int r = 0; r < L.outerSize(); ++r){

            double a_ii = L.coeff(r, r);

            // Skip boundary value (FIXME: we should rely on "std::vector boundaries" )
            if( a_ii == 1.0 )
                continue;

            double sum = 0.0;
            double mean = 0.0;
            // Look up columns
            //std::cout << "r: " << r << std::endl;
            for(Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(L, r); it; ++it) {
                if( it.col() == r)
                    continue;

                //std::cout << "value: " << it.value() << std::endl;
                sum += harmonic_weight_map[it.col()] * it.value();
                mean += it.value();
            }

            // Right hand side is always null in our case
            sum =  (0/*rhs[r]*/-sum) / (a_ii);

            buffer[r] = harmonic_weight_map[r]*(1.0 - t) + sum*t;
        }
        buffer.swap( harmonic_weight_map);
    }

    if( (iter)%2 != 0 ){
        std::cout << "swap" << std::endl;
        buffer.swap( harmonic_weight_map );
    }

    return;
}

// -----------------------------------------------------------------------------
