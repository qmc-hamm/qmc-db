export function fmtErr(value, err) {
  if (typeof value !== "number" || !Number.isFinite(value)) return "—";
  const dec2 = () => value.toFixed(Math.abs(value) >= 1e21 ? 0 : 2);
  if (typeof err !== "number" || !Number.isFinite(err) || err <= 0) return dec2();
  let d = Math.floor(Math.log10(err));
  if (-d > 20) return dec2();
  let digit = Math.round(err / 10 ** d);
  if (digit >= 10) { d += 1; digit = Math.round(err / 10 ** d); }
  if (d > 0) {
    const q = 10 ** d;
    return (Math.round(value / q) * q).toFixed(0) + "(" + digit * q + ")";
  }
  return value.toFixed(-d) + "(" + digit + ")";
}

export function fmtET(isoUtc) {
  const s = new Date(isoUtc).toLocaleString("en-US", {
    timeZone: "America/New_York", year: "numeric", month: "short", day: "numeric",
    hour: "numeric", minute: "2-digit", hour12: true,
  });
  return s + " ET";
}

// Multi-select filter with an "all" option: true → keep only "all" (clear the rest).
export function keepAll(wasAll, allNow, othersNow) {
  return (allNow && !wasAll) || !othersNow;
}
