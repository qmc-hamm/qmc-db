import test from "node:test";
import assert from "node:assert/strict";
import { fmtErr, fmtET, keepAll } from "../static/fmt.js";

test("fmtErr", () => {
  assert.equal(fmtErr(-85.144, 0.016), "-85.14(2)");
  assert.equal(fmtErr(0.08559, 0.00012), "0.0856(1)");
  assert.equal(fmtErr(-1430.7127, null), "-1430.71");
  assert.equal(fmtErr(NaN, null), "—");
  assert.equal(fmtErr(null, 0.1), "—");
  assert.equal(fmtErr(12.3, 0.96), "12(1)");
  assert.equal(fmtErr(1234, 30), "1230(30)");
  assert.equal(fmtErr(1.5, NaN), "1.50");
});

test("fmtET", () => {
  const s = fmtET("2026-11-19T18:00:00Z");
  assert.ok(s.includes("1:00"), s);
  assert.ok(s.includes("PM"), s);
  assert.ok(s.endsWith("ET"), s);
});

test("fmtErr guards", () => {
  assert.equal(fmtErr(1.5, 0), "1.50");
  assert.equal(fmtErr(1.5, Infinity), "1.50");
  assert.equal(fmtErr(1.5, -0.1), "1.50");
  assert.equal(fmtErr("abc", 0.1), "—");
  assert.equal(fmtErr(Infinity, 0.1), "—");
  assert.equal(fmtErr(1.5, 1e-30), "1.50");
});

test("keepAll", () => {
  assert.equal(keepAll(false, true, true), true);   // clicked "all" → clears others
  assert.equal(keepAll(true, true, true), false);   // added a value → drops "all"
  assert.equal(keepAll(false, false, false), true); // nothing left → back to "all"
  assert.equal(keepAll(false, false, true), false);
});
