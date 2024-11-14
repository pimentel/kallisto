#ifndef BLOCK_ARRAY_H
#define BLOCK_ARRAY_H

#include <vector>
#include <algorithm>

template<typename T = int>
struct block {
    uint32_t lb, ub;
    T val;

    block() : lb(0), ub(0) {}
    block(uint32_t idx) : lb(idx) {}
    block(uint32_t lb, uint32_t ub, T val) : lb(lb), ub(ub), val(val) {}
    ~block() {}

    bool operator<(const block<T>& rhs) const {
        return (this->lb < rhs.lb);
    }
};

template<typename T = int>
class BlockArray {
    private:

    public:
        union { // TODO: could optimize memory if just vector w/o union
            block<T> mono;
            std::vector<block<T> > poly;
        };
        uint8_t flag;

        BlockArray() : flag(0), mono(0) {
        }

        ~BlockArray() {
            switch(flag) {
                case 0:
                case 1:
                    mono.~block();
                    break;
                case 2:
                    poly.~vector();
            }
        }

        BlockArray(const BlockArray& o) {

            flag = o.flag;
            switch(flag) {
                case 0:
                case 1:
                    new (&mono) auto(o.mono);
                    break;
                case 2:
                    new (&poly) auto(o.poly);
            }
        }

        BlockArray(BlockArray&& o) {

            flag = o.flag;
            o.flag = 0;

            switch(flag) {
                case 0:
                case 1:
                    new (&mono) auto(o.mono);
                    o.mono = block<T>();
                    break;
                case 2:
                    new (&poly) auto(o.poly);
                    o.poly.clear();
            }
        }

        BlockArray& operator=(const BlockArray& o) {

            clear();
            flag = o.flag;
            switch(flag) {
                case 0:
                case 1:
                    new (&mono) auto(o.mono);
                    break;
                case 2:
                    new (&poly) auto(o.poly);
            }
            return *this;
        }

        BlockArray& operator=(BlockArray&& o){

            if (this != &o) {

                flag = o.flag;
                o.flag = 0;

                switch(flag) {
                    case 0:
                    case 1:
                        new (&mono) auto(o.mono);
                        o.mono = block<T>();
                        break;
                    case 2:
                        new (&poly) auto(o.poly);
                        o.poly.clear();
                }
            }

            return *this;
        }

        void insert(uint32_t lb, uint32_t ub, T val) {

            // TODO: Make sure block does not overlap with existing block
            if (flag == 0) {
                mono.lb = lb;
                mono.ub = ub;
                mono.val = val;
                flag = 1;

                return;
            } else if (flag == 1) {
                block<T> m = mono;
                mono.~block();
                new (&poly) std::vector<block<T> >;
                poly.push_back(m);
                flag = 2;
            }

            poly.emplace_back(lb, ub, val);
            std::sort(poly.begin(), poly.end());
        }

        void insert(block<T>& b) {

            if (flag == 0) {
                mono = b;
                flag = 1;
                return;
            } else if (flag == 1) {
                block<T> m = mono;
                mono.~block();
                new (&poly) std::vector<block<T> >;
                poly.push_back(m);
                flag = 2;
            }
            poly.push_back(b);
            std::sort(poly.begin(), poly.end());

        }

        void overwrite(uint32_t lb, uint32_t ub, T& val) {
            // Overwrites part of the BlockArray with a new value

            if (flag == 1) {
                // If the unitig is monochrome

                if (mono.lb < lb) {
                    // lb is an internal index

                    block<T> ow(lb, ub, val);
                    block<T> m = mono;
                    if (mono.ub == ub) {
                        // [ mono ][ overwrite ]
                        m.ub = lb;
                        mono.~block();
                        new (&poly) std::vector<block<T> >;
                        // Pushing blocks in-order: don't need to sort
                        poly.push_back(m);
                        poly.push_back(ow);
                        flag = 2;
                    } else {
                        // [ mono ][ overwrite ][ mono ]

                        block<T> rest(ub, m.ub, m.val);
                        m.ub = lb;
                        mono.~block();
                        new (&poly) std::vector<block<T> >;
                        // Pushing blocks in-order: don't need to sort
                        poly.push_back(m);
                        poly.push_back(ow);
                        poly.push_back(rest);
                        flag = 2;
                    }

                } else if (mono.ub > ub) {
                    // ub is an internal index
                    // [ overwrite ][ mono ]

                    block<T> ow(lb, ub, val);
                    block<T> m = mono;
                    m.lb = ub;
                    mono.~block();
                    new (&poly) std::vector<block<T> >;
                    poly.push_back(m);
                    poly.push_back(ow);
                    flag = 2;
                } else {
                    // [ overwrite ]
                    // mono.lb == lb, mono.ub == ub
                    mono.val = val;

                }

            } else {
                // If the unitig is polychrome, we need to find all the blocks
                // that overlap [lb, ub]

                if (poly[0].lb == lb && poly[poly.size()-1].ub == ub) {
                    // [ overwrite ]
                    // Overwrite spans entire unitig
                    poly.~vector();

                    new (&mono) block<T>(lb, ub, val);
                    flag = 1;

                    return;

                }

                std::vector<block<T> > rep;
                for (auto& b : poly) {
                    if (b.ub < lb || b.lb > ub) {
                        // ... [ b ] ... [ overwrite ] ...
                        // or
                        // ... [ overwrite ] ... [ b ] ...
                        rep.push_back(b);
                    } else if (b.lb >= lb && b.ub <= ub) {
                        // ... [ overwrite [ b ] overwrite ] ...
                        continue;
                    } else if (b.lb < lb) {
                        // ... [ b ] overwrite ] ...

                        rep.emplace_back(b.lb, lb, b.val);
                        if (b.ub > ub) {
                            // ... [ b [ overwrite ] b ] ...
                            rep.emplace_back(ub, b.ub, b.val);
                        }
                    } else if (b.ub > ub) {
                        // ... [ overwrite [ b ] ...
                        rep.emplace_back(ub, b.ub, b.val);
                    }
                }

                // Finally insert new block
                rep.emplace_back(lb, ub, val);
                std::sort(rep.begin(), rep.end());
                poly.clear();
                poly.assign(rep.begin(), rep.end());

            }

        }

        const T& operator[](size_t idx) const {

            if (flag < 2) {
                return mono.val;
            }

            // TODO: Find upper bound without having to create temp struct
            block<T> tmp(idx);
            auto ub = std::upper_bound(poly.begin(), poly.end(), tmp);
            // if (poly.blocks.size() > 0) {
            //     block<T> bb = poly.blocks[poly.blocks.size() - 1];
            // }

            if (ub == poly.begin()) {
                // No elements with lb <= idx
                // TODO: Handle this gracefully
                std::cout << "not found????" << std::endl;
            }
            return (--ub)->val;
        }

        std::pair<uint32_t, T> block_index(uint32_t idx) const {
            if (flag < 2) {
                return std::make_pair(0, mono.val);
            }

            block<T> tmp(idx);
            auto ub = std::upper_bound(poly.begin(), poly.end(), tmp);
            --ub;
            return std::make_pair(ub - poly.begin(), (ub)->val);
        }

        // Returns the ECs of the block up to and including the EC at idx
        std::vector<T> get_leading_vals(size_t idx) const {
            std::vector<T> vals;
            if (flag == 1) {
                vals.push_back(mono.val);
            } else if (flag == 2) {
                vals.reserve(poly.size());
                for (const auto& b : poly) {
                    if (b.lb <= idx) {
                        vals.push_back(b.val);
                    }
                }
            }

            return vals;
        }

        // Returns the ECs of the block up to and including the EC at idx
        std::vector<T> get_trailing_vals(size_t idx) const {
            std::vector<T> vals;
            if (flag == 1 && mono.lb >= idx) {
                vals.push_back(mono.val);
            } else if (flag == 2) {
                vals.reserve(poly.size());
                for (const auto& b : poly) {
                    if (b.lb >= idx) {
                        vals.push_back(b.val);
                    }
                }
            }

            return vals;
        }

        std::pair<uint32_t, uint32_t> get_block_at(size_t idx) const {
            if (flag == 1) {
                return std::make_pair(mono.lb, mono.ub);
            }

            block<T> tmp(idx);
            auto ub = std::upper_bound(poly.begin(), poly.end(), tmp);

            if (ub == poly.begin()) {
                // No elements with lb <= idx
                return std::make_pair(-1, -1);
            }

            --ub;

            return std::make_pair(ub->lb, ub->ub);
        }

        void reserve(uint32_t sz) {
            if (sz < 2) return;
            flag = 2;
            poly.reserve(sz);
        }

        void clear() {
            switch(flag) {
                case 0:
                case 1:
                    mono.~block();
                    break;
                case 2:
                    poly.~vector();
            }
            flag = 0;
        }

        size_t size() const {
            if (flag < 2) return flag;
            return poly.size();
        }

        size_t length() const {
            if (flag < 2) return mono.ub;
            return poly[poly.size()-1].ub;
        }

        void get_vals(std::vector<T>& vals) const {
            vals.clear();
            if (flag == 1) {
                vals.push_back(mono.val);
            } else if (flag == 2){
                vals.reserve(poly.size());
                for (const auto& b : poly) {
                    vals.push_back(b.val);
                }
            }
        }

        typename std::vector<block<T> >::iterator begin() {
            return poly.begin();
        }

        typename std::vector<block<T> >::iterator end() {
            return poly.end();
        }

        typename std::vector<block<T> >::const_iterator begin() const {
            return poly.begin();
        }

        typename std::vector<block<T> >::const_iterator end() const  {
            return poly.end();
        }

        void print() const {

            if (flag == 0) return;
            if (flag == 1) {
                std::cout << "[" << mono.lb << ", " << mono.ub << "): " << mono.val;
            } else {
                for (const auto& b : poly) {
                    std::cout << "[" << b.lb << ", " << b.ub << "): " << b.val << ", ";
                }
            }
            std::cout << std::endl;
        }

        // Extracts the slice [lb, ub) from *this and shifts it such that it
        // occupies the range [0, ub-lb).
        BlockArray<T> get_slice(uint32_t lb, uint32_t ub) const {
            BlockArray<T> slice;

            if (flag == 1) {
                slice.insert(lb, ub, mono.val);
            } else if (flag == 2) {

                block<T> tmp(lb);
                auto upper = std::upper_bound(poly.begin(), poly.end(), tmp);

                if (upper == poly.begin()) {
                    // No elements with lower bound <= lb
                    // TODO: Handle this
                }

                --upper;

                do {
                    slice.insert(
                            std::max(lb, upper->lb)-lb,
                            std::min(ub, upper->ub)-lb,
                            upper->val
                    );
                    ++upper;
                } while(upper != poly.end() && upper->lb < ub);
            }

            return slice;
        }


        void serialize(std::ostream& out) const {

            out.write((char *)&flag, sizeof(flag));
            if (flag == 0) return;
            else if (flag == 1) {

                out.write((char *)&mono.lb, sizeof(mono.lb));
                out.write((char *)&mono.ub, sizeof(mono.ub));
                mono.val.serialize(out);
            } else {

                size_t tmp_size = poly.size();
                out.write((char *)&tmp_size, sizeof(tmp_size));
                for (const auto b : poly) {
                    out.write((char *)&b.lb, sizeof(b.lb));
                    out.write((char *)&b.ub, sizeof(b.ub));
                    b.val.serialize(out);
                }
            }
        }

        void deserialize(std::istream& in, bool additional_opt) {

            clear();

            in.read((char *)&flag, sizeof(flag));
            if (flag == 0) {
                return;
            } else if (flag == 1) {

                in.read((char *)&mono.lb, sizeof(mono.lb));
                in.read((char *)&mono.ub, sizeof(mono.ub));
                T val;
                val.deserialize(in, additional_opt);
                mono.val = std::move(val);
            } else {

                size_t tmp_size;
                uint32_t lb, ub;
                in.read((char *)&tmp_size, sizeof(tmp_size));
                mono.~block();
                new (&poly) std::vector<block<T> >;
                poly.reserve(tmp_size);

                for (size_t i = 0; i < tmp_size; ++i) {
                    in.read((char *)&lb, sizeof(lb));
                    in.read((char *)&ub, sizeof(ub));

                    T val;
                    val.deserialize(in, additional_opt);
                    insert(lb, ub, std::move(val));
                }
            }
        }
};

#endif // BLOCK_ARRAY_H
