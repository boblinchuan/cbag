// SPDX-License-Identifier: BSD-3-Clause AND Apache-2.0
/*
Copyright (c) 2018, Regents of the University of California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Copyright 2019 Blue Cheetah Analog Design Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef CBAG_COMMON_VIA_PARAM_H
#define CBAG_COMMON_VIA_PARAM_H

#include <array>

#include <cbag/common/typedefs.h>
#include <cbag/common/vector.h>

namespace cbag {
namespace layout {

struct via_param {
  public:
    std::array<cnt_t, 2> num = {0, 0};
    vector cut_dim = {0, 0};
    vector cut_spacing = {0, 0};
    std::array<vector, 2> enc = {vector{0, 0}, vector{0, 0}};
    std::array<vector, 2> off = {vector{0, 0}, vector{0, 0}};
    int priority_ = 1;

    via_param();

    via_param(cnt_t vnx, cnt_t vny, offset_t w, offset_t h, offset_t vspx, offset_t vspy,
              offset_t enc1l, offset_t enc1r, offset_t enc1t, offset_t enc1b, offset_t enc2l,
              offset_t enc2r, offset_t enc2t, offset_t enc2b, int priority);

    bool operator==(const via_param &rhs) const noexcept;

    std::string to_string() const noexcept;
};

} // namespace layout
} // namespace cbag

#endif
