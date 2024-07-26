#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop.hxx>
#include <loop_device.hxx>

#include <AMReX.H>

#include <setup_eos.hxx>

#include "linear_interp_ND.hxx"

#include <vector>

namespace TestMe {

using namespace amrex;
using namespace std;
using namespace Loop;

extern "C" void Tester(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTSX_Tester;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("EOS Tester");

  std::vector<double> input_x{};
  std::vector<double> input_y{};
  std::vector<double> control_x{};
  std::vector<double> control_y{};

  double start{-5.5};
  double end{5.5};
  double step{0.2};
  double step_control{0.001};

  for ( int i = 0; i < static_cast<int>( (end-start)/step  ); ++i ) {
  	input_x.push_back(i*step+start);
          input_y.push_back( pow(i*step+start,2) );
  }
  
  for ( int i = 0; i < static_cast<int>( (end-start)/step_control  ); ++i ) {
  	control_x.push_back(i*step_control+start);
          control_y.push_back( pow(i*step_control+start,2) );
  }
  
  const size_t n = input_x.size();
  
  auto x_ptr = std::unique_ptr<double[]>(new double[n]);
  auto y_ptr = std::unique_ptr<double[]>(new double[n]);
  
  for (int i = 0; i < n; ++i) x_ptr[i] = std::move(input_x[i]);
  for (int i = 0; i < n; ++i) y_ptr[i] = std::move(input_y[i]);
  
  std::array<size_t,1> num_points = {n};
  
  linear_interp_uniform_ND_t<double, 1, 1>  interpolator = linear_interp_uniform_ND_t<double, 1, 1>(
        std::move(y_ptr), std::move(num_points), std::move(x_ptr));

  auto *ptr = &interpolator;

  cctk_grid.loop_all_device<1, 1, 1>(
        cctk_grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

	CCTK_REAL xx = p.x;
        dumdum(p.I) = ptr->interpolate<0>(xx)[0];
	//dumdum(p.I) = interpolator.interpolate<0>(xx)[0];

	});

}

} // namespace TestMe
