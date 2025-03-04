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

#include <fmt/core.h>

#include <cbag/layout/via_param.h>

namespace cbag {
namespace layout {

void set_via_enc_offset(vector &enc, vector &off, offset_t encl, offset_t encr, offset_t enct,
                        offset_t encb) {
    enc[0] = (encr + encl) / 2;
    enc[1] = (enct + encb) / 2;
    off[0] = encr - enc[0];
    off[1] = enct - enc[1];
}

via_param::via_param() = default;

via_param::via_param(cnt_t vnx, cnt_t vny, offset_t w, offset_t h, offset_t vspx, offset_t vspy,
                     offset_t enc1l, offset_t enc1r, offset_t enc1t, offset_t enc1b, offset_t enc2l,
                     offset_t enc2r, offset_t enc2t, offset_t enc2b, int priority)
    : num{vnx, vny}, cut_dim{{w, h}}, cut_spacing{vspx, vspy} {
    set_via_enc_offset(enc[0], off[0], enc1l, enc1r, enc1t, enc1b);
    set_via_enc_offset(enc[1], off[1], enc2l, enc2r, enc2t, enc2b);
    priority_ = priority;
}

bool via_param::operator==(const via_param &rhs) const noexcept {
    return num == rhs.num && cut_dim == rhs.cut_dim && cut_spacing == rhs.cut_spacing &&
           enc == rhs.enc && off == rhs.off && priority_ == rhs.priority_;
}

std::string via_param::to_string() const noexcept {
    return fmt::format(
        "ViaParam(num=({}, {}), dim=({}, {}), sp=({}, {}), enc=[({}, {}), ({}, {})], "
        "off=[({}, {}), ({}, {}))], priority_={})",
        num[0], num[1], cut_dim[0], cut_dim[1], cut_spacing[0], cut_spacing[1], enc[0][0],
        enc[0][1], enc[1][0], enc[1][1], off[0][0], off[0][1], off[1][0], off[1][1], priority_);
}

} // namespace layout
} // namespace cbag
