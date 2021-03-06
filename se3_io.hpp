// Copyright 2011 Ethan Eade. All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:

//    1. Redistributions of source code must retain the above
//       copyright notice, this list of conditions and the following
//       disclaimer.

//    2. Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials
//       provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY ETHAN EADE ``AS IS'' AND ANY EXPRESS
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ETHAN EADE OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

// The views and conclusions contained in the software and
// documentation are those of the authors and should not be
// interpreted as representing official policies, either expressed or
// implied, of Ethan Eade.

#ifndef LATL_SE3_IO_HPP
#define LATL_SE3_IO_HPP

#include <latl/se3.hpp>
#include <latl/io.hpp>

namespace latl
{
    template <class Scalar>
    std::ostream& operator<<(std::ostream& out, const SE3<Scalar>& se3)
    {
        for (int i=0; i<3; ++i)
            out << se3.rotation().matrix()[i]
                << Vector<1,Scalar>(se3.translation()[i]) << std::endl;
        return out;
    }

    template <class Scalar>
    std::istream& operator>>(std::istream& in, SE3<Scalar>& se3)
    {
        Matrix<3,4,Scalar> Rt;
        if (in >> Rt)
            se3 = SE3<Scalar>(Rt.T()[3], SO3<Scalar>::from_matrix(slice<0,0,3,3>(Rt)));
        return in;
    }    
}

#endif
