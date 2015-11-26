/*  Copyright 2014 IST Austria

Contributed by: Jan Reininghaus, Stefan Huber

This file is part of DIPHA.

DIPHA is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DIPHA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with DIPHA.  If not, see <http://www.gnu.org/licenses/>. */

#include <dipha/includes.h>
#ifdef USE_FGT
#include "figtree.h"
#endif

void print_help_and_exit()
{
    std::cerr << "Usage: " << "diagram_distance [options] --dim N --time T input_filename_1 input_filename_2 ..." << std::endl;
    std::cerr << std::endl;
    std::cerr << "--distance_squared    --  print square of norm of difference" << std::endl;
    std::cerr << "--inner_product       --  print inner product" << std::endl;
    std::cerr << "--dim N               --  target dimension" << std::endl;
    std::cerr << "--time T              --  heat diffusion time" << std::endl;
    std::cerr << "--help                --  prints this screen" << std::endl;
#ifdef USE_FGT
    std::cerr << "--use_fgt             --  use Fast-Gauss-Transform" << std::endl;
#endif
    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
}

#ifdef USE_FGT
void parse_command_line( int argc, char** argv, bool& distance_squared, bool& inner_product, bool& use_fgt, int64_t& dim, double& time, std::vector< std::string >& input_filenames )
#else
void parse_command_line( int argc, char** argv, bool& distance_squared, bool& inner_product, int64_t& dim, double& time, std::vector< std::string >& input_filenames )
#endif
{
#ifdef USE_FGT
    use_fgt = false;
#endif
    distance_squared = false;
    inner_product = false;
    
    if( argc < 6 )
        print_help_and_exit();

    for( int idx = 1; idx < argc; idx++ ) {
        const std::string option = argv[ idx ];
        if( option == "--help" ) {
            print_help_and_exit( );
        } else if( option == "--distance_squared" ) {
            distance_squared = true;
        } else if( option == "--inner_product" ) {
            inner_product = true;
#ifdef USE_FGT
        } else if (option == "--use_fgt") {
            use_fgt = true;
#endif
        } else if( option == "--dim" ) {
            idx++;
            if( idx >= argc )
                print_help_and_exit( );
            std::string parameter = std::string( argv[ idx ] );
            size_t pos_last_digit;
            dim = std::stoll( parameter, &pos_last_digit );
            if( pos_last_digit != parameter.size( ) )
                print_help_and_exit( );
        } else if( option == "--time" ) {
            idx++;
            if( idx >= argc )
                print_help_and_exit( );
            std::string parameter = std::string( argv[ idx ] );
            size_t pos_last_digit;
            time = std::stod( parameter, &pos_last_digit );
            if( pos_last_digit != parameter.size( ) )
                print_help_and_exit( );
        } else {
            input_filenames.push_back( argv[ idx ] );
        }
    }
}

#ifdef USE_FGT
double L2_diagram_inner_product_fgt( const std::vector< double >& births_1, const std::vector< double >& deaths_1,
                                 const std::vector< double >& births_2, const std::vector< double >& deaths_2, double time )
{
    assert( births_1.size( ) == deaths_1.size( ) );
    assert( births_2.size( ) == deaths_2.size( ) );

    std::vector< double > mu1;
    std::vector< double > mu2;

    int64_t num_points_1 = (int64_t)births_1.size( );
    for( int64_t idx = 0; idx < num_points_1; idx++ ) {
        mu1.push_back( births_1[ idx ] );
        mu1.push_back( deaths_1[ idx ] ); 
    }

    int64_t num_points_2 = (int64_t)births_2.size( );
    for( int64_t idx = 0; idx < num_points_2; idx++ ) {
        mu2.push_back( births_2[ idx ] );
        mu2.push_back( deaths_2[ idx ] );    
    }

    double *g_auto = new double[ num_points_2 ];    // IGFT output
    double local_sum = 0.0;                         // Inner product
    double bw = std::sqrt(2.0*time);                // Bandwidth
                                  
    int64_t local_begin = dipha::element_distribution::get_local_begin( births_1.size() );
    int64_t local_end = dipha::element_distribution::get_local_end( births_1.size() );

    std::vector< double > source;
    for( int64_t idx = local_begin; idx < local_end; idx++ ) {
        source.push_back( mu1[ 2*idx+0 ] );
        source.push_back( mu1[ 2*idx+1 ] );
    }

    std::vector< double > weight( source.size()/2, 1.0 ); // Weights   

    // Call to IFGT
    memset( g_auto, 0, sizeof( double )*births_2.size( ) );
    figtree( 
        2,                // Dimensionality
        source.size()/2,  // Nr. of points in source  
        num_points_2,     // Nr. of points in target
        1,                // Nr. of weights 
        &source[0],       // Source points
        bw,               // Bandwidth
        &weight[0],       // Weights
        &mu2[0],          // Target points
        1e-5,             // Tolerance 
        g_auto );

    for( int i = 0; i < births_2.size(); i++) {
        local_sum += g_auto[i];
    }

    // Mirror source points by swapping the coordinates in source
    for( int64_t idx = 0; idx < source.size()/2; idx++ ) {
        std::swap( source[ 2*idx+0 ], source[ 2*idx+1 ] );
    }

    // Call to IGFT for mirrored points
    memset( g_auto, 0, sizeof( double )*births_2.size( ) );
    figtree( 
        2,                // Dimensionality
        source.size()/2,  // Nr. of points in source  
        num_points_2,     // Nr. of points in target
        1,                // Nr. of weights 
        &source[0],       // Source points
        bw,               // Bandwidth
        &weight[0],       // Weights
        &mu2[0],          // Target points
        1e-5,             // Tolerance 
        g_auto );

    for( int i = 0; i < births_2.size(); i++) {
        local_sum -= g_auto[i];
    }

    // Gather and sum
    std::vector< double > gathered_values;
    dipha::mpi_utils::all_gather( local_sum, gathered_values );
    double sum = std::accumulate( gathered_values.begin( ), gathered_values.end( ), 0.0, std::plus< double >( ) );
    
    // Cleanup
    delete [] g_auto;

    // Return final kernel value
    const double pi = 3.14159265358979323846;
    return (sum / ( 2.0 * pi * time ));
}
#endif

double L2_diagram_inner_product( const std::vector< double >& births_1, const std::vector< double >& deaths_1,
                                 const std::vector< double >& births_2, const std::vector< double >& deaths_2, double time )
{
    std::vector< double > a;
    std::vector< double > mu_1;
    std::vector< double > mu_2;

    std::vector< double > b;
    std::vector< double > nu_1;
    std::vector< double > nu_2;

    // Sanity check
    assert( births_1.size( ) == deaths_1.size( ) );
    assert( births_2.size( ) == deaths_2.size( ) );

    int64_t num_points_1 = (int64_t)births_1.size( );
    for( int64_t idx = 0; idx < num_points_1; idx++ ) {
        a.push_back( 1.0 );
        mu_1.push_back( births_1[ idx ] );
        mu_2.push_back( deaths_1[ idx ] );

        // mirrored points
        a.push_back( -1.0 );
        mu_1.push_back( deaths_1[ idx ] );
        mu_2.push_back( births_1[ idx ] );
    }

    int64_t num_points_2 = (int64_t)births_2.size( );
    for( int64_t idx = 0; idx < num_points_2; idx++ ) {
        b.push_back( 1.0 );
        nu_1.push_back( births_2[ idx ] );
        nu_2.push_back( deaths_2[ idx ] );

        // mirrored points
        b.push_back( -1.0 );
        nu_1.push_back( deaths_2[ idx ] );
        nu_2.push_back( births_2[ idx ] );
    }


    double local_sum = 0.0;
    int64_t local_begin = dipha::element_distribution::get_local_begin( 2 * num_points_1 );
    int64_t local_end = dipha::element_distribution::get_local_end( 2 * num_points_1 );
    for( int64_t idx_1 = local_begin; idx_1 < local_end; idx_1++ ) {
        for( int64_t idx_2 = 0; idx_2 < 2 * num_points_2; idx_2++ ) {
            double mu_nu_squared = ( mu_1[ idx_1 ] - nu_1[ idx_2 ] ) * ( mu_1[ idx_1 ] - nu_1[ idx_2 ] )
                + ( mu_2[ idx_1 ] - nu_2[ idx_2 ] ) * ( mu_2[ idx_1 ] - nu_2[ idx_2 ] );
            local_sum += a[ idx_1 ] * b[ idx_2 ] * std::exp( -mu_nu_squared / ( 2 * time ) );
        }
    }

    std::vector< double > gathered_values;
    dipha::mpi_utils::all_gather( local_sum, gathered_values );
    double sum = std::accumulate( gathered_values.begin( ), gathered_values.end( ), 0.0, std::plus< double >( ) );
    
    const double pi = 3.14159265358979323846;
    double sum_R2 = sum / ( 2.0 * pi * time );
    double half_plane = sum_R2 / 2.0;
    return half_plane;
}

double L2_diagram_distance_squared( const std::vector< double >& births_1, const std::vector< double >& deaths_1,
                                const std::vector< double >& births_2, const std::vector< double >& deaths_2, double time )
{
    std::vector< double > births = births_1;
    std::vector< double > deaths =  deaths_1;

    // HACK: mirrored diagram results in subtraction
    std::copy( deaths_2.begin( ), deaths_2.end(), std::back_inserter( births ) );
    std::copy( births_2.begin( ), births_2.end( ), std::back_inserter( deaths ) );

    return L2_diagram_inner_product( births, deaths, births, deaths, time );
}


void read_diagram( const std::string& filename, std::vector< int64_t >& dims, std::vector< double >& births, std::vector< double >& deaths )
{
    MPI_File file = dipha::mpi_utils::file_open_read_only( filename );

    // read preamble
    std::vector< int64_t > preamble;
    dipha::mpi_utils::file_read_at_vector( file, 0, 3, preamble );

    int64_t num_points = preamble[ 2 ];

    dims.resize( num_points );
    births.resize( num_points );
    deaths.resize( num_points );

    std::vector< int64_t > int64_t_temp;
    dipha::mpi_utils::file_read_at_vector( file, preamble.size( ) * sizeof( int64_t ), num_points * 3, int64_t_temp );

    std::vector< double > double_temp;
    dipha::mpi_utils::file_read_at_vector( file, preamble.size( ) * sizeof( int64_t ), num_points * 3, double_temp );

    for( int64_t cur_point = 0; cur_point < num_points; cur_point++ ) {
        dims[ cur_point ] = int64_t_temp[ 3 * cur_point + 0 ];
        births[ cur_point ] = double_temp[ 3 * cur_point + 1 ];
        deaths[ cur_point ] = double_temp[ 3 * cur_point + 2 ];
    }

    MPI_File_close( &file );
}

void get_selected_points( const std::string& filename, int64_t dim, std::vector< double >& births, std::vector< double >& deaths )
{
    std::vector< double > all_births;
    std::vector< double > all_deaths;
    std::vector< int64_t > all_dims;
    read_diagram( filename, all_dims, all_births, all_deaths );

    for( int64_t idx = 0; idx < (int64_t)all_dims.size( ); idx++ ) {
        if( all_dims[ idx ] == dim ) {
            births.push_back( all_births[ idx ] );
            deaths.push_back( all_deaths[ idx ] );
        }
    }
}

int main( int argc, char** argv )
{
    // mandatory MPI initilization call
    MPI_Init( &argc, &argv );

    std::vector< std::string > arg_filenames; 
    int64_t dim;
    double time;
    bool distance_squared;
    bool inner_product;
#ifdef USE_FGT
    bool use_fgt; 
    parse_command_line( argc, argv, distance_squared, inner_product, use_fgt, dim, time, arg_filenames );
#else
    parse_command_line( argc, argv, distance_squared, inner_product, dim, time, arg_filenames );
#endif    
    std::vector< std::string > input_filenames;
    for( const auto& filename : arg_filenames ) {
        if( !dipha::file_types::is_dipha_file( filename  ) ) {
            std::ifstream input_stream( filename.c_str() );
            if( input_stream.fail() ) {
                dipha::mpi_utils::error_printer_if_root() << filename << " cannot be read!" << std::endl;
                MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
            }
            std::string line;
            while( !getline( input_stream, line ).eof() ) {
                input_filenames.push_back( line );
                line.clear();
            }
        } else if( dipha::file_types::get_file_type( filename ) == dipha::file_types::PERSISTENCE_DIAGRAM ) {
            input_filenames.push_back( filename );
        } else {
            dipha::mpi_utils::error_printer_if_root( ) << filename << " is a DIPHA file but not a PERSISTENCE_DIAGRAM!" << std::endl;
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
        }
    }

    int64_t num_datasets = (int64_t) input_filenames.size( );
    std::vector< std::vector< double > > gram_matrix( num_datasets );
    for( auto& row : gram_matrix )
        row.resize( num_datasets );
    
    for( int64_t idx_1 = 0; idx_1 < num_datasets; idx_1++ ) {
        const std::string& input_filename_1 = input_filenames[ idx_1 ];
        std::vector< double > births_1;
        std::vector< double > deaths_1;
        get_selected_points( input_filename_1, dim, births_1, deaths_1 );

        // loop over upper triangular part (without diagonal)
        for( int64_t idx_2 = idx_1 + 1; idx_2 < num_datasets; idx_2++ ) {
            //std::cout << "idx_1:" << idx_1 << ", idx_2: " << idx_2 << std::endl;
            const std::string& input_filename_2 = input_filenames[ idx_2 ];
            std::vector< double > births_2;
            std::vector< double > deaths_2;
            get_selected_points( input_filename_2, dim, births_2, deaths_2 );
   
            double value = 0.0;
            if( inner_product )
#ifdef USE_FGT
                if ( use_fgt )
                    value = L2_diagram_inner_product_fgt( births_1, deaths_1, births_2, deaths_2, time );
                else
                    value = L2_diagram_inner_product( births_1, deaths_1, births_2, deaths_2, time );                    
#else
                value = L2_diagram_inner_product( births_1, deaths_1, births_2, deaths_2, time );
#endif            
            else {
                double squared = L2_diagram_distance_squared( births_1, deaths_1, births_2, deaths_2, time );
                if( distance_squared )
                    value = squared;
                else
                    value = std::sqrt( squared );
            }
            gram_matrix[ idx_1 ][ idx_2 ] = value;
            gram_matrix[ idx_2 ][ idx_1 ] = value;
        }
        
        // fill diagonal
        if( inner_product )
            gram_matrix[ idx_1 ][ idx_1 ] = L2_diagram_inner_product( births_1, deaths_1, births_1, deaths_1, time );
        else
            gram_matrix[ idx_1 ][ idx_1 ] = 0.0;
    }

    // print matrix to cout
    dipha::mpi_utils::cout_if_root( ) << std::scientific << std::setprecision( 16 );
    for( const auto& row : gram_matrix ) {
        for( const auto& value : row )
            dipha::mpi_utils::cout_if_root() << value << " ";
        dipha::mpi_utils::cout_if_root( ) << std::endl;
    }
   
    //// REGRESSION TEST:
    //std::vector< double > births_1;
    //std::vector< double > deaths_1;
    //births_1.push_back( 0 );
    //deaths_1.push_back( 1 );
    //std::vector< double > births_2;
    //std::vector< double > deaths_2;
    //births_2.push_back( 2 );
    //deaths_2.push_back( 3 );
    //time = 1;
    //dipha::mpi_utils::cout_if_root( ) << L2_diagram_distance_squared( births_1, deaths_1, births_2, deaths_2, time );
    ////==> distance_squared = 0.197525 

    MPI_Finalize();
}