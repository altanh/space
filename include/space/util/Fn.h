#pragma once

#include <array>
#include <tuple>
#include <utility>

namespace space {

namespace fn {

/// List library

template <size_t... _Is> using List = std::index_sequence<_Is...>;
using Nil = std::index_sequence<>;
template <size_t _N> using CountTo = std::make_index_sequence<_N>;

template <size_t _I, typename _Seq> struct Tail;
template <size_t _Head, size_t... _Tail> struct Tail<0, List<_Head, _Tail...>> {
  using T = List<_Head, _Tail...>;
};

template <> struct Tail<0, Nil> {
  using T = Nil;
};

template <size_t _I, size_t _Head, size_t... _Tail>
struct Tail<_I, List<_Head, _Tail...>> {
  using T = typename Tail<_I - 1, List<_Tail...>>::T;
};

template <typename _S1, typename _S2> struct Concat;

template <size_t... _S1, size_t... _S2>
struct Concat<List<_S1...>, List<_S2...>> {
  using T = List<_S1..., _S2...>;
};

template <size_t _Start, size_t _End> struct Range;
template <size_t _End> struct Range<_End, _End> {
  using T = Nil;
};

template <size_t _Start, size_t _End> struct Range {
  using T =
      typename Concat<List<_Start>, typename Range<_Start + 1, _End>::T>::T;
};

static_assert(std::is_same_v<typename Range<0, 3>::T, CountTo<3>>);

template <size_t _N, typename _Seq> struct Take;

template <typename _Seq> struct Take<0, _Seq> {
  using T = Nil;
};

template <size_t _Head, size_t... _Tail> struct Take<1, List<_Head, _Tail...>> {
  using T = List<_Head>;
};

template <size_t _N, size_t _Head, size_t... _Tail>
struct Take<_N, List<_Head, _Tail...>> {
  using T =
      typename Concat<List<_Head>, typename Take<_N - 1, List<_Tail...>>::T>::T;
};

template <size_t _N, typename _Seq> struct Contains;

template <size_t _N> struct Contains<_N, Nil> {
  static constexpr bool v = false;
};

template <size_t _N, size_t _Head, size_t... _Tail>
struct Contains<_N, List<_Head, _Tail...>> {
  static constexpr bool v = _N == _Head || Contains<_N, List<_Tail...>>::v;
};

/// Tuple manipulation library

// https://stackoverflow.com/a/58674921
template <class F, class T, size_t N = 0>
void ApplyAtDyn(const T &t, size_t i, const F &f) {
  if (N == i) {
    f(std::integral_constant<size_t, N>(), std::get<N>(t));
    return;
  }
  if constexpr (N + 1 < std::tuple_size_v<T>) {
    return ApplyAtDyn<F, T, N + 1>(t, i, f);
  }
}

template <typename F, typename... Ts, size_t... Indices>
decltype(auto) Apply(const std::tuple<Ts...> &t, F &&f, List<Indices...>) {
  // if f(t[i]) is non-void, return nothing
  if constexpr ((std::is_void_v<decltype(std::forward<F>(f)(
                     std::get<Indices>(t)))> ||
                 ...)) {
    (std::forward<F>(f)(std::get<Indices>(t)), ...);
  } else {
    return std::make_tuple(std::forward<F>(f)(std::get<Indices>(t))...);
  }
}

template <typename F, typename... Ts>
decltype(auto) Apply(const std::tuple<Ts...> &t, F &&f) {
  return Apply(t, std::forward<F>(f), CountTo<sizeof...(Ts)>{});
}

template <typename F, typename... Ts, size_t... Indices>
decltype(auto) ApplyEnumerate(const std::tuple<Ts...> &t, F &&f,
                              List<Indices...>) {
  // if f(t[i]) is non-void, return nothing
  if constexpr ((std::is_void_v<decltype(std::forward<F>(f)(
                     std::integral_constant<size_t, Indices>(),
                     std::get<Indices>(t)))> ||
                 ...)) {
    (std::forward<F>(f)(std::integral_constant<size_t, Indices>(),
                        std::get<Indices>(t)),
     ...);
  } else {
    return std::make_tuple(std::forward<F>(f)(
        std::integral_constant<size_t, Indices>(), std::get<Indices>(t))...);
  }
}

template <typename F, typename... Ts>
decltype(auto) ApplyEnumerate(const std::tuple<Ts...> &t, F &&f) {
  return ApplyEnumerate(t, std::forward<F>(f), CountTo<sizeof...(Ts)>{});
}

template <typename F, typename... Ts, size_t... Indices>
bool ApplyUntil(const std::tuple<Ts...> &t, F &&f, List<Indices...>) {
  return (std::forward<F>(f)(std::get<Indices>(t)) || ...);
}

template <typename F, typename... Ts, size_t... Indices>
bool ApplyEnumerateUntil(const std::tuple<Ts...> &t, F &&f, List<Indices...>) {
  return (std::forward<F>(f)(std::integral_constant<size_t, Indices>(),
                             std::get<Indices>(t)) ||
          ...);
}

template <typename F, typename... Ts>
bool ApplyEnumerateUntil(const std::tuple<Ts...> &t, F &&f) {
  return ApplyEnumerateUntil(t, std::forward<F>(f), CountTo<sizeof...(Ts)>{});
}

template <typename F, typename... Ts>
bool ApplyUntil(const std::tuple<Ts...> &t, F &&f) {
  return ApplyUntil(t, std::forward<F>(f), CountTo<sizeof...(Ts)>{});
}

template <typename... Ts, size_t... Is>
decltype(auto) Select(const std::tuple<Ts...> &t, List<Is...>) {
  return std::make_tuple(std::get<Is>(t)...);
}

/// Array read/write library

template <typename V, size_t Out, size_t In, size_t... Dst, size_t... Src>
void Copy(std::array<V, In> &out, const std::array<V, Out> &in, List<Dst...>,
          List<Src...>) {
  (..., (out[Dst] = in[Src]));
}

// Generic version?
template <typename TOut, typename TIn, size_t... Dst, size_t... Src>
void Copy(TOut &out, const TIn &in, List<Dst...>, List<Src...>) {
  (..., (std::get<Dst>(out) = std::get<Src>(in)));
}

template <typename V, size_t Out, size_t In, size_t... Dst>
void Scatter(std::array<V, Out> &out, const std::array<V, In> &in,
             List<Dst...> dst) {
  static_assert(sizeof...(Dst) == In);
  Copy(out, in, dst, CountTo<In>{});
}

template <typename V, size_t Out, size_t In, size_t... Src>
void Gather(std::array<V, Out> &out, const std::array<V, In> &in,
            List<Src...> src) {
  static_assert(sizeof...(Src) == Out);
  Copy(out, in, CountTo<Out>{}, src);
}

} // namespace fn

} // namespace space