// SPDX-License-Identifier: BSD-3-Clause
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
*/

#ifndef CBAG_NETLIST_NGSPICE_H
#define CBAG_NETLIST_NGSPICE_H

#include <string>
#include <unordered_map>

#include <cbag/netlist/core.h>
#include <cbag/netlist/nstream_output.h>

namespace cbag {
namespace netlist {

class ngspice_stream : public nstream_output {
  private:
    cnt_t precision_ = 6;
    std::unordered_map<std::string, std::string> cell_name_map_;
    std::unordered_map<std::string, std::string> prop_name_map_;
    std::unordered_map<std::string, std::unordered_map<std::string, std::string>> cell_prop_map_;

  public:
    ngspice_stream();

    ngspice_stream(const std::string &fname, cnt_t precision);

    const std::string &get_cell_name(const std::string &cell_name) const;

    const std::string &get_prop_name(const std::string &cell_name,
                                     const std::string &prop_name) const;

    cnt_t precision() const noexcept;
};

template <> struct traits::nstream<ngspice_stream> {
    using type = ngspice_stream;

    static void close(type &stream);

    static void write_header(type &stream, const std::vector<std::string> &inc_list, bool shell);

    static void write_end(type &stream);

    static void write_cv_header(type &stream, const std::string &name,
                                const sch::cellview_info &info, bool shell, bool write_subckt,
                                bool write_declarations);

    static void write_cv_end(type &stream, const std::string &name, bool write_subckt);

    static void write_unit_instance(type &stream, const std::string &prefix, cnt_t inst_idx,
                                    const spirit::ast::name_unit &name_ast,
                                    const term_net_vec_t &conn_list, const param_map &params,
                                    const sch::cellview_info &info,
                                    const net_rename_map_t *net_map_ptr);

    static void append_netlist(type &stream, const std::string &netlist);

    static void write_supply_wrapper(type &stream, const std::string &name,
                                     const sch::cellview_info &info);

    static net_rename_map_t new_net_rename_map();
};

} // namespace netlist
} // namespace cbag

#endif // CBAG_NETLIST_ngspice_H
