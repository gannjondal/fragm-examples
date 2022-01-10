#info 3D Newton implementation for the polynomial q^n + JuliaC=0 with real JuliaC

//See the original thread at https://fractalforums.org/fractal-mathematics-and-new-theories/28/revisiting-the-3d-newton/1026
//Basic implementation by gannjondal

/** Background, and explanation:
  - Meaning of the identifiers in below explanations:   
    + z stands for the variable that is visible to the fractal program, and that is used for bailout check, calculation of the distance estimation etc etc   
    + q stands for the variable that is used for the calculation of the Newton method   
   
  - Explanation of the "trick".   
    It was necessary to keep the following conditions in mind:   
    + First of all: To get a useful 3D image it is necessary to define an 'inside', and an 'outside' region (keep in mind that a clean Newton fractal converges everywhere)
    + This flavor of Newton 3D fractal implementation has originally been implemented for MB3D. 
	  The common bailout control, and DE mechanisms of MB3D, (and also MB2) need a DIVERGENT iteration value (here z)    
    + It is of course necessary to have a variable (here q) that is used for the actual Newton iteration
    + In MB3D (and maybe also in MB2) it is not possible to keep another triplex value than z between the iterations.   
   
    Therefore it was necessary to have a z that diverges - AND it must be possible to always calculate q from z uniquely, and vice versa 
	(mathematically: There must be a one-to-one relation between q, and z).   
    Furthermore I decided that the area that converges to one specific value (below 'solution') is taken as 'outside'. The other areas are 'inside'.     
    Of course that way you see only the outer surface of the 'inside' area, but it turned out that with some tweaks it's possible to see ... more of the transition area.
   
    Below algorithm fullfills all above conditions - Each iteration looks like:   
    a)  Extract q[n] as q[n] = 1/z[n] + solution   
    b)  Run the Newton method calculation which provdes a new q[n+1]   
    c)  Calculate z[n+1] as z[n+1] = 1/(q[n+1]-solution), and hand it over to the fractal program   
   
    To have a classical Newton star one would need to skip the very first inversion.   
    But I decided to keep that first inversion -   
    It transforms the formula to a finite object - so to speak a Newton bulb.   
    And this object is of course much more handleable than the classical infinite star - in fact I see that as a benefit which has been introduced without any cost....   
**/

/** Parameters to manipulate the Newton calculation:
Power:   
   Default = 3   
   The power of the polynom   
    (x,y,z)^power - c = 0   
   to be solved   
  
Fac_Phi, Fac_Theta:   
   Default = 1 for both   
   With this params you can vary the radial angles of the potence when calculating the power of (x,y,z).   
   Hence it disturbs the power calculation, but in a very smooth way.   
   You even may try to set one of the params to 0.   
	
Vary_Inversion:
   Default = -1
   On a first view that looks similar like the Fake_Bailout_xyz params.
   But it's for y, and z axes only - and they are not equivalent to x here.
   Try values not too far from -1 (start with -0.75 ... -1.5)
  
Fake_Bailout_pre, Fake_Bailout_post:   
   Default = 1 / 1   
   ORIGINAL INTENTION:   
     In some hybrids it is necessary to increase the value of the bailout on the formula window to avoid unwanted cuts.   
     In this case you may set both Fake_Bailout_xyz to bailout/4 as a start value.   
	 It is especially amazing to use different values for both - try for instance Fake_Bailout_pre = 1, Fake_Bailout_post = -1
   MEANWHILE:   
     It turned out that changing this value can lead to completely different structures.   
     I don't have a simple suggestion - you may use samples that exist especially on ff.org   
   ATTENTION:   
     Changing the value of both params in parallel scales the complete object.   
     You may add a general scaling (when used in hybrids) to scale the object back.   
   Fake_Bailout_pre - used BEFORE the transformation   
   Fake_Bailout_post - used AFTER the transformation   
**/

#include "MathUtils.frag"
#include "DE-Raytracer.frag"

#group NewtonVarPower
   //Power n for the equation q^n + JuliaC = 0
   uniform float Power; slider[-32,3,32]
   
   //Known solution of the polynomial (and origin of the implicite inversion)
   uniform vec3 Solution; slider[(-4,-4,-4),(1,0,0),(4,4,4)]
   //HINT:  Start with Solution = 1.0, Julia=true, JuliaC=(-1,0,0)
   //       Other params for Solution, and JuliaC are mainly useful if using this frag as pretransform
   
   //Core parameters to tweak the fractal
   //Factor to vary the phi angle of the power calculations
   uniform float Factor_Phi; slider[-20,1,20]
   //Factor to vary the theta angle of the power calculations
   uniform float Factor_Theta; slider[-20,1,20]
   //Vary the calculation of the absolute value of q BEFORE the Newton calc (HINT:  Try also negative values) 
   uniform float Fake_Bailout_pre; slider[-20,1,20]
   //Vary the calculation of the absolute value of q AFTER the Newton calc (HINT:  Try also negative values) 
   uniform float Fake_Bailout_post; slider[-20,1,20]
   //Vary the inversion of q (similar to Fake_Bailout, but affects only y, and z)
   uniform float Vary_Inversion; slider[-4,-1,0]
   
   //Core iteration controls - FOR NEWTON + FOLDING ONLY
   //Offset to avoid the division by zero, you usually don't need to change the value
   uniform float Offset; slider[1.0e-15,0.0001,1]
   //Usual bailout value (for |1/(q-solution)|). ATTENTION:  Needs to be changed for powers <> 3!
   uniform float Bailout; slider[0.001,4,1024]
   //Number of Newton iterations. HINT: Start with a value of 1 or 2 - here the Newton function is used as pretransform !!
   uniform int Iterations; slider[0,1,200]
   //Number of Newton iteration taken into account for color calculation
   uniform int ColorIterations; slider[0,3,200]
   
   //DE tweaks - Scale
   uniform float DEscale1; slider[0,1,2]
   //DE tweaks - Offset
   uniform float DEoffset1; slider[0,0,8]
   
   //Julia + c controls
   //HINT:  Start with Solution = 1.0, Julia=true, JuliaC=(-1,0,0)
   //       Other params for Solution, and JuliaC are mainly useful if using this frag as pretransform
   //Julia calculation?
   uniform bool Julia; checkbox[true]
   //C value for Julia calculation
   uniform vec3 JuliaC; slider[(-4,-4,-4),(-1,0,0),(4,4,4)]

float sq_r, r, r1, dr, theta, phi, r_pow, theta_pow, phi_pow, pow_eff, fac_eff, cth, cph, sph, sth, tmpx, tmpy, tmpz, tmpx2, tmpy2, tmpz2, r_zxy, r_cxy, h, scale;
vec3 c;
int i;
  
float DE(vec3 pos) {

// Preparation operations
    vec3 z = pos;
	r1 = length(z);
	
    pow_eff = 1.0 - Power;
    fac_eff = (Power - 1.0)/Power;
    
    dr = 1.0;
    scale = 1.0;
    i = 0;
	
    c = (Julia ? JuliaC : pos);
	r_cxy = sqrt(c.x*c.x + c.y*c.y);

//Main iteration loop
    while(r1<Bailout && (i<Iterations)) {
       dr = dr*r1*2.0;
	   
    // Converting the diverging (x,y,z) back to the variable
    // that can be used for the (converging) Newton method calculation
          sq_r = Fake_Bailout_pre/(dot(z,z) + Offset);
          z.x = z.x*sq_r + Solution.x;
          z.y = Vary_Inversion*z.y*sq_r + Solution.y;
          z.z = Vary_Inversion*z.z*sq_r + Solution.z;
           
    // Calculate 1/z^(power-1)
          sq_r = dot(z,z);
          r = sqrt(sq_r);
          phi = Factor_Phi*asin(z.z/r) ;
          theta = Factor_Theta*atan(z.y,z.x) ;
            
          r_pow = pow(r, pow_eff);
          phi_pow = phi*pow_eff;
          theta_pow = theta*pow_eff;
      
          cth = cos(theta_pow);
          sth = sin(theta_pow);
          cph = cos(phi_pow);
          sph = sin(phi_pow);
      
          tmpx = -cph*cth*r_pow/Power;
          tmpy = -cph*sth*r_pow/Power;
          tmpz = sph*r_pow/Power;
      
    // Multiply c and z
          r_zxy = sqrt(tmpx*tmpx + tmpy*tmpy);
          
          h = 1 - c.z*tmpz/(r_cxy*r_zxy + Offset);
          
          tmpx2 = (c.x*tmpx - c.y*tmpy)*h;
          tmpy2 = (c.y*tmpx + c.x*tmpy)*h;
          tmpz2 = r_cxy*tmpz + r_zxy*c.z;
   
    // Bring everything together        
          z.x = fac_eff*z.x + tmpx2;
          z.y = fac_eff*z.y + tmpy2;
          z.z = fac_eff*z.z + tmpz2;
         
    // Below the hack that provides a divergent value of (x,y,z) to MB3D
    // although the plain Newton method does always converge
          sq_r = Fake_Bailout_post/((dot(z-Solution,z-Solution))+ Offset);
          
          z.x = (z.x - Solution.x)*sq_r;
          z.y = Vary_Inversion*(z.y - Solution.y)*sq_r;
          z.z = Vary_Inversion*(z.z - Solution.z)*sq_r;
           
    //DE helper calculations (?)
          dr = dr * DEscale1 + DEoffset1;
    //     dr =  pow(r, pow_eff)*(Power - 1.0)*dr + 10.0;
    //     dr = max(dr*DerivativeBias, r_pow*dr*pow_eff + 1.0);
          r1 = length(z);
        
          if (i<ColorIterations) orbitTrap = min(orbitTrap, abs(vec4(z.x,z.y,z.z,r1*r1)));
          i++;
    }
   return 0.5*log(r1)*r1/dr;
}

#preset CommonNewtonBulb
FOV = 0.4
Eye = 1.41821122,1.39130151,-2.070901
Target = -5.39926173,-2.73930773,3.9673178
Up = -0.09550437,0.868544433,0.486322093
EquiRectangular = false
AutoFocus = false
FocalPlane = 30
Aperture = 0
Gamma = 2.08335
ToneMapping = 3
Exposure = 0.6522
Brightness = 1
Contrast = 1
AvgLumin = 0.5,0.5,0.5
Saturation = 1
LumCoeff = 0.212500006,0.715399981,0.0720999986
Hue = 0
GaussianWeight = 1
AntiAliasScale = 2
DepthToAlpha = false
ShowDepth = false
DepthMagnitude = 1
Detail = -4.5
DetailAO = -5
FudgeFactor = 0.75
MaxDistance = 3000
MaxRaySteps = 2500
Dither = 0
NormalBackStep = 0 NotLocked
AO = 0.203921571,0.203921571,0.203921571,0.83
Specular = 0.6
SpecularExp = 14
SpecularMax = 20
SpotLight = 1,0.90196079,0.745098054,1
SpotLightDir = 0.6,0.5
CamLight = 1,1,1,1.5
CamLightMin = 0.07
Glow = 1,1,1,0.5
GlowMax = 52
Fog = 0
HardShadow = 0.35 NotLocked
ShadowSoft = 12.6
QualityShadows = true
Reflection = 0 NotLocked
DebugSun = false NotLocked
BaseColor = 0.760784328,0.760784328,0.760784328
OrbitStrength = 0.4
X = 1,1,1,1
Y = 0.345097989,0.666666985,0,0.02912
Z = 1,0.666666985,0,1
R = 0.0784313977,1,0.941175997,-0.0194
BackgroundColor = 0.713725507,0.866666675,0.596078455
GradientBackground = 0.3
CycleColors = true
Cycles = 5
EnableFloor = false NotLocked
FloorNormal = 0,0,0
FloorHeight = 0
FloorColor = 1,1,1
Power = 3
Solution = 1,0,0
Factor_Phi = 1
Factor_Theta = 1
Fake_Bailout_pre = 1
Fake_Bailout_post = 1
Vary_Inversion = -1
Offset = 1e-5
Bailout = 4
Iterations = 40
ColorIterations = 7
DEscale1 = 1
DEoffset1 = 4
Julia = true
JuliaC = -1,0,0
#endpreset

#preset NewtonTree
FOV = 0.4
Eye = 1.64284857,1.57141761,-2.34849137
Target = -3.78684741,-2.70384153,4.87926882
Up = -0.078253214,0.882715894,0.463345206
EquiRectangular = false
AutoFocus = false
FocalPlane = 30
Aperture = 0
Gamma = 2.08335
ToneMapping = 3
Exposure = 0.6522
Brightness = 1
Contrast = 1
AvgLumin = 0.5,0.5,0.5
Saturation = 1
LumCoeff = 0.212500006,0.715399981,0.0720999986
Hue = 0
GaussianWeight = 1
AntiAliasScale = 2
DepthToAlpha = false
ShowDepth = false
DepthMagnitude = 1
Detail = -5.25
DetailAO = -5
FudgeFactor = 1
MaxDistance = 100
MaxRaySteps = 2500
Dither = 0.5
NormalBackStep = 5 NotLocked
AO = 0.203921571,0.203921571,0.203921571,0.83
Specular = 0.6
SpecularExp = 14
SpecularMax = 20
SpotLight = 1,0.90196079,0.745098054,1
SpotLightDir = 0.6,0.5
CamLight = 1,1,1,1.5
CamLightMin = 0.07
Glow = 1,1,1,0.5
GlowMax = 52
Fog = 0
HardShadow = 0.35 NotLocked
ShadowSoft = 12.6
QualityShadows = true
Reflection = 0 NotLocked
DebugSun = false NotLocked
BaseColor = 0.760784328,0.760784328,0.760784328
OrbitStrength = 0.4
X = 1,1,1,1
Y = 0.345097989,0.666666985,0,0.02912
Z = 1,0.666666985,0,1
R = 0.0784313977,1,0.941175997,-0.0194
BackgroundColor = 0.713725507,0.866666675,0.596078455
GradientBackground = 0.3
CycleColors = true
Cycles = 5
EnableFloor = false NotLocked
FloorNormal = 0,0,0
FloorHeight = 0
FloorColor = 1,1,1
Power = 4
Solution = 1,0,0
Factor_Phi = 1
Factor_Theta = 1
Fake_Bailout_pre = 1
Fake_Bailout_post = -1
Vary_Inversion = -1
Offset = 1e-5
Bailout = 11
Iterations = 40
ColorIterations = 7
DEscale1 = 1
DEoffset1 = 4
Julia = true
JuliaC = -1,0,0
#endpreset

#preset NewtonBulbColorVariant
FOV = 0.9472675
Eye = 0.590536351,0.393237328,-1.14650676
Target = -0.191921449,0.19073232,-0.541399037
Up = -0.28661842,-0.119444438,-0.410596174
EquiRectangular = false
AutoFocus = false
FocalPlane = 7
Aperture = 0
Gamma = 0.4863813
ToneMapping = 5
Exposure = 1
Brightness = 1
Contrast = 1
AvgLumin = 0.5,0.5,0.5
Saturation = 1
LumCoeff = 0.212500006,0.715399981,0.0720999986
Hue = 0
GaussianWeight = 1
AntiAliasScale = 1
DepthToAlpha = false
ShowDepth = false
DepthMagnitude = 1
Detail = -4.5
DetailAO = -0.19
FudgeFactor = 0.2
MaxDistance = 100
MaxRaySteps = 4821
Dither = 0
NormalBackStep = 0 NotLocked
AO = 0.203921571,0.203921571,0.203921571,0.8
Specular = 0.05
SpecularExp = 20
SpecularMax = 10
SpotLight = 1,0.937254902,0.858823529,0.5
SpotLightDir = 0.80952382,0.52380954
CamLight = 1,0.90196079,0.831372559,0.95
CamLightMin = 0.075
Glow = 1,1,1,0.5
GlowMax = 25
Fog = 0
HardShadow = 0.85 NotLocked
ShadowSoft = 0
QualityShadows = true
Reflection = 0.337254912 NotLocked
DebugSun = false NotLocked
BaseColor = 1,1,1
OrbitStrength = 0.3
X = 0.207843095,0.600000024,0,1
Y = 0.00392156886,0.235294104,1,1
Z = 0.996078372,1,0.725490212,1
R = 1,0.0862745121,0.00392156886,1
BackgroundColor = 0.168627456,0.368627459,0.474509805
GradientBackground = 1.5
CycleColors = true
Cycles = 2.75
EnableFloor = true NotLocked
FloorNormal = 0,0,1
FloorHeight = 1.569485
FloorColor = 0.0784313679,0.176470593,0.333333313
Power = 3
Solution = 1,0,0
Factor_Phi = 1
Factor_Theta = 1
Fake_Bailout_pre = 1
Fake_Bailout_post = 1
Vary_Inversion = -1
Offset = 1e-5
Bailout = 4
Iterations = 40
ColorIterations = 7
DEscale1 = 1
DEoffset1 = 4
Julia = true
JuliaC = -1,0,0
#endpreset