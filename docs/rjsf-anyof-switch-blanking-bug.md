# RJSF anyOf field-blanking bug: switching away and back clears values

**Applies to:** any workflow using an `anyOf`-based "Method"-style selector (a
Union of pydantic models rendered via RJSF's `AnyOfField`/`MultiSchemaField`),
where the schema is authored via pydantic + `wt_compiler`/`rjsf-overrides` and
rendered by ecoscope-desktop's RJSF (`@rjsf/utils` v5.19.4).

**Case study:** `home_range_args.args` in this repo (`spec.yaml`), a
`McpMethodArgs | EtdMethodArgs | BbmmMethodArgs` selector
(`wd-partner-tasks/.../tasks/_home_range.py`, `_bbmm.py`).

---

## 1. The symptom

Selecting a Method (e.g. ETD) and filling in its fields works fine. Switching
to a different Method (e.g. MCP) and back to the first one **blanks every
field that was previously populated** - not reset to defaults, genuinely
empty. The bug compounds: the first visit to any option is always fine, but
every option becomes permanently blank-on-revisit the moment you switch away
from it once.

One easy way to misread this bug: if a field happens to be **object-typed**
(e.g. a nested `anyOf` sub-selector), it can *look* immune because its first
option coincides with its own default (so the "reset" is invisible) - this
is a rendering coincidence, not real protection, and does not generalize.

## 2. Root cause (confirmed via source, not guessed)

Do not trust behavioral guesses about RJSF here - trace the actual source.
The relevant files (`@rjsf/utils` v5.19.4, fetch from
`https://raw.githubusercontent.com/rjsf-team/react-jsonschema-form/v5.19.4/packages/utils/src/...`):

- `schema/sanitizeDataForNewSchema.ts`
- `schema/getDefaultFormState.ts` (`computeDefaults`, `maybeAddDefaultToObject`)
- `mergeDefaultsWithFormData.ts`
- `schema/retrieveSchema.ts`
- `schema/getClosestMatchingOption.ts`, `getMatchingOption.ts`, `getFirstMatchingOption.ts`
- `packages/core/src/components/fields/MultiSchemaField.tsx` (this is the
  actual `AnyOfField` implementation)

When the user picks a new option, `MultiSchemaField.onOptionChange` runs:

```js
let newFormData = schemaUtils.sanitizeDataForNewSchema(newOption, oldOption, formData);
if (newFormData && newOption) {
  newFormData = schemaUtils.getDefaultFormState(newOption, newFormData, 'excludeObjectChildren');
}
onChange(newFormData, undefined, this.getFieldId());
```

**`sanitizeDataForNewSchema`** does NOT delete a field that no longer applies
- it sets it to `undefined` as an explicit **own-property** (a tombstone),
not a removed key:

```js
Object.keys(oldSchema.properties).forEach((key) => {
  if (has(data, key)) removeOldSchemaData[key] = undefined;
});
```

That tombstone rides forward in `formData` on every subsequent switch,
because the final assembly is a spread: `{ ...data, ...removeOldSchemaData, ...nestedData }`.

**`mergeDefaultsWithFormData`** (called from `getDefaultFormState`) then
does this for objects:

```js
const acc = Object.assign({}, defaults);           // starts correct
return Object.keys(formData).reduce((acc, key) => {
  acc[key] = mergeDefaultsWithFormData(get(defaults, key), get(formData, key));
  return acc;
}, acc);
```

and its scalar/leaf case is:

```js
return formData;   // unconditional - discards a correctly computed default
```

Since `Object.keys(formData)` includes tombstoned keys (present, just
`undefined`-valued), any field that has ever been switched away from gets its
freshly computed default **stomped back to `undefined` on every future visit
to that option.** This is the entire bug.

### Two non-buggy code paths this abuses to fix itself

`sanitizeDataForNewSchema` has two paths that do NOT tombstone:

1. **Same name + same type, both old and new schema define it** → direct
   value swap, no tombstone:
   ```js
   if (newOptionDefault !== formValue) {
     if (oldOptionDefault === formValue) removeOldSchemaData[key] = newOptionDefault;
     // else: value is left untouched (user's customization is preserved)
   }
   ```
2. **Object-typed key with no old-schema counterpart** → recurses; if the
   recursive call returns `undefined` (nothing to carry forward), the key is
   **left absent** rather than added as a tombstone:
   ```js
   if (newSchemaTypeForKey === 'object' || ...) {
     const itemData = sanitizeDataForNewSchema(newKeyedSchema, oldKeyedSchema, formValue);
     if (itemData !== undefined || newSchemaTypeForKey === 'array') {
       nestedData[key] = itemData;   // only added if itemData is real
     }
   }
   ```
   Critically, this recursive call only returns something non-`undefined`
   (even an empty `{}`) if the schema being recursed into has a
   **`properties` key of its own**. A bare `anyOf` (no sibling `properties`)
   returns `undefined` and gets the same non-immunity as a scalar.

Branch **matching** (which option is "selected") is a separate mechanism
(`MultiSchemaField` constructor / `componentDidUpdate` →
`schemaUtils.getClosestMatchingOption` → `getMatchingOption` first, falling
back to `calculateIndexScore` scoring by `default`/`const` match). This does
**not** require disjoint property names across options - it only needs each
option to have at least one field that scores distinctly for its own data
(typically its own required disambiguator). Confirmed by hand-computing the
score for this repo's exact schema.

## 3. The fix - two techniques, chosen per field

Both route every field through one of the two non-buggy paths above instead
of through the tombstone-prone "bare scalar with no counterpart" path.

### Technique A - shared property names

For any field meaning the same thing across ≥2 variants (here: `crs`,
`percentiles`, `expansion_factor` for ETD/BBMM), give it the **same literal
property name and type** across all variants that share the concept. This
triggers the direct-value-swap repair - no tombstone ever created for it.

Caveat: this only protects a switch **directly** between two schemas that
both define the name. If a third variant that lacks the field sits in
between (e.g. `Etd → Mcp → Etd`, and Mcp has no `expansion_factor`), the
tombstone reappears. If more than 2 of your N variants share a field, that
field still needs Technique B for full protection (see `expansion_factor`
below).

### Technique B - wrap genuinely-unique fields in a tiny object

For a field unique to one variant (or shared by fewer than all variants),
wrap it in its own single-field pydantic model, e.g.:

```python
class EtdSpeedSettings(BaseModel):
    model_config = ConfigDict(extra="ignore")
    max_speed_factor: Annotated[float, Field(title="Max Speed Factor (Kilometers per Hour)", ...)]

class EtdMethodArgs(BaseModel):
    speed_settings: Annotated[EtdSpeedSettings, Field(title="Max Speed Factor (Kilometers per Hour)", json_schema_extra={...})]
    ...
```

This gives the field real `properties` of its own, so once it's tombstoned
it self-heals (absent, not tombstoned, on the next switch away-and-back).

**Do not use a docstring on the wrapper class** - pydantic renders a model's
docstring as its schema `description`, which leaks any internal
implementation notes into the rendered form as if they were user-facing
help text. Use a plain `#` comment above the class instead.

### The part that actually took two rounds to get right: `allOf` vs bare `$ref`

Wrapping alone is not sufficient. Any `Annotated[Model, Field(title=..., ...)]`
- i.e. a model-typed field with a title override or *any* other `Field(...)`
metadata alongside it - makes pydantic emit:

```json
{"allOf": [{"$ref": "#/$defs/EtdSpeedSettings"}], "default": {...}, "title": "..."}
```

(JSON Schema's standard workaround for "no siblings alongside `$ref`".)

But `sanitizeDataForNewSchema` only resolves a **bare** `{"$ref": ...}`
(`has(schema, '$ref')`, checked directly - never one nested inside `allOf`).
So the wrapper field's own `type` is **never** resolved to `'object'`, and it
silently falls through to the untyped-scalar (tombstone-prone) path -
defeating the wrapping entirely, with no error or warning anywhere. This
affects every wrapped field here, including genuinely-required ones - a
`title=` override alone is enough to trigger the `allOf` wrap regardless of
required-ness.

Two ways to fix this were tried and rejected before landing on the one that
works - **both failures were confirmed by direct reproduction**, not assumed:

1. **Inject `type`/`properties` into the field's own `json_schema_extra`,
   computed from `model_cls.model_json_schema()`.** This actually works for
   a field with a real default and no internal `$ref`. But:
   - If the wrapped model's own inner field is itself `$ref`-bearing (e.g.
     wrapping a nested `anyOf`), pydantic's schema generator raises a
     `KeyError` on **any** raw `"$ref"` string placed inside
     `json_schema_extra` - it walks the whole generated schema for `$ref`
     bookkeeping and can't reconcile one it didn't discover via its own
     type-driven generation.
   - Separately, pydantic **silently drops** an injected `"type"` key
     specifically when the referencing field has **no real python default**
     (i.e. is genuinely required) - every other injected key survives, only
     `type` does not. This blocks exactly the fields that need to stay
     required for backend disambiguation.
2. Because of (1), this cannot be fixed purely in Python for every field.

**What actually works:** inject `type`/`properties` via `spec.yaml`'s
`rjsf-overrides.$defs` section instead - the same override mechanism this
file already uses elsewhere (e.g. `ValueGrouper.properties.index_name.oneOf`,
`Group Data.groupers.groupers.items.anyOf`). This operates on the **final
assembled schema**, after pydantic's own generation is already done, where
neither restriction applies:

```yaml
rjsf-overrides:
  $defs:
    EtdMethodArgs.properties.speed_settings.type: object
    EtdMethodArgs.properties.speed_settings.properties:
      max_speed_factor:
        default: 1.05
        description: "..."
        title: Max Speed Factor (Kilometers per Hour)
        type: number
```

The `properties` value here is a **hand-copied duplicate** of the wrapper
model's own real schema (get the exact values by inspecting the compiled
`rjsf.json`'s `$defs.<WrapperModel>.properties` after a first
`regenerate_rjsf.sh` pass) - not auto-derived, since the auto-derivation
approach is exactly what breaks in Python. Keep the Python wrapper model as
the source of truth and treat this override as a manually-synced mirror;
update both together if the wrapped field ever changes.

### One more field-shape gotcha: title placement and font size

RJSF renders an **object's own title** as a big fieldset-style header,
distinct from a scalar field's normal-sized `<label>`. Once a field is
wrapped in a model (Technique B), its outer wrapper is object-typed and its
title would render oversized and visually inconsistent with the rest of the
form. Fix: hide the **wrapper's** own label via `uiSchema`
(`ui:options.label: false`) and put the real, user-facing title on the
**inner scalar field** instead - its normal scalar rendering keeps the same
font size as every other (non-wrapped) field on the form. This is the
opposite of what seems natural at first (hiding the inner field and keeping
the outer title, since the outer title is what already existed pre-wrapping)
- get this backwards and every wrapped field renders with a jarring
oversized heading.

Note this does **not** apply to `type: array` fields (e.g. a multi-select
"Percentile Levels") - RJSF renders arrays with the same big heading style
as objects regardless of wrapping, since both are "container" types as
opposed to true scalars. That is pre-existing RJSF/ecoscope-desktop
behavior, unrelated to this fix, and not something fixable from schema
authoring alone (it lives in ecoscope-desktop's `ArrayFieldTemplate`).

## 4. Step-by-step recipe for a new `anyOf` selector

1. **Identify fields that mean the same thing across ≥2 variants.** Give
   them the exact same property name and JSON type across every variant
   that shares the concept (Technique A). If fewer than *all* variants share
   it, also apply Technique B to it (see `expansion_factor` in this repo -
   shared by ETD/BBMM but not MCP, so it's *also* wrapped on both sides with
   matching wrapper shape).
2. **Wrap every field still unique to one variant** in its own tiny
   single-field pydantic model (Technique B). No docstring on the wrapper
   class - use a `#` comment instead.
3. Give each wrapper field a `title=` matching what the field's label should
   say to the user - this is safe regardless of the `allOf` issue in step 4,
   since the fix is applied downstream in `spec.yaml`.
4. Run `./dev/regenerate_rjsf.sh` once, then inspect the compiled
   `rjsf.json`'s `$defs.<WrapperModel>.properties` for the exact shape to
   hand-copy.
5. Add one `$defs` override pair per wrapped field in `spec.yaml`:
   ```yaml
   <ParentModel>.properties.<wrapper_field>.type: object
   <ParentModel>.properties.<wrapper_field>.properties:
     <inner_field_name>: <hand-copied from step 4>
   ```
6. Add a `uiSchema` override per wrapped field: hide the wrapper's own
   label, leave the inner field's label visible.
   ```yaml
   uiSchema:
     Home Range Method.home_range_args.args.anyOf:
       - <option 0>
         <wrapper_field>:
           ui:options:
             label: false
   ```
7. `./dev/regenerate_rjsf.sh` again and verify via step 5 below.
8. Each variant keeps (or gets) at least one genuinely-required field (real
   or wrapped) with no python default, so branch matching/backend
   disambiguation for a bare `param.yaml` submission still works. The one
   variant with *no* required field at all is the "fallback" - verify only
   one variant fits that description.

## 5. How to verify a fix like this (don't trust it without doing this)

Do **not** trust a fix based on reasoning about RJSF's behavior alone, and do
**not** trust a simulation built by hand-flattening schemas before feeding
them in - both of those produced false confidence during this
investigation, twice. The only reliable verification method:

1. Fetch the real `@rjsf/utils` source for the exact version pinned in
   `ecoscope-desktop`'s `package.json` (`v5.19.4` here) - do not paraphrase
   it, port it verbatim.
2. Port `sanitizeDataForNewSchema`, `computeDefaults`/`getDefaultFormState`,
   `mergeDefaultsWithFormData` faithfully, preserving quirks like the
   `includeUndefinedValues === 'excludeObjectChildren'` special value.
3. Extract the **actual compiled `rjsf.json`** (via
   `./dev/regenerate_rjsf.sh`) - not a hand-written test schema - and
   recursively resolve its `$ref`/single-element-`allOf` the same way RJSF's
   own `retrieveSchema` would (this recursive resolution is what caught the
   `allOf` vs bare-`$ref` distinction; skipping it is what produced the
   first false-positive "fix").
4. Simulate `MultiSchemaField.onOptionChange`'s exact call sequence -
   `sanitizeDataForNewSchema` then `getDefaultFormState(..., 'excludeObjectChildren')`
   - across a full cycle through **every** variant, **more than once each**
   (a 2-variant back-and-forth is not enough - the `expansion_factor`
   MCP-in-the-middle bug only showed up on a 3-variant loop). After each
   step, assert every property the new schema defines is not `undefined`.
5. Only after the harness passes, hand the resulting standalone schema +
   uiSchema fragment to a human to spot-check in
   `https://app.ecoscope.io/configuration-form-playground` (or the real
   app) - this is what caught the first "fix" not actually working live,
   despite the (at-the-time incomplete) harness reporting success.

The harness script and driver used for this investigation are not committed
to this repo (they lived in a scratch directory) - reconstruct them from
this document if needed; the port is straightforward once you have the
exact RJSF source files listed in §2.

## 6. Files touched by the actual fix in this repo

- `wd-partner-tasks/src/ecoscope-workflows-ext-wd/ecoscope_workflows_ext_wd/tasks/_home_range.py`
  - `EtdMethodArgs`, `McpMethodArgs`: shared `crs`/`percentiles`.
  - `EtdSpeedSettings`, `EtdGridCellSizeSettings`, `EtdExpansionFactorSettings`,
    `McpRingsCorrectionSettings`: wrapper models.
  - Module-level comment documents the full history (superseded disjoint-name
    design kept inline, marked `SUPERSEDED`, for context).
- `wd-partner-tasks/src/ecoscope-workflows-ext-wd/ecoscope_workflows_ext_wd/tasks/_bbmm.py`
  - `BbmmMethodArgs`: shared `crs`/`percentiles`/`expansion_factor`.
  - `BbmmLocationErrorSettings`, `BbmmTimeStepSettings`,
    `BbmmMaxDataGapSettings`, `BbmmExpansionFactorSettings`: wrapper models.
- `spec.yaml`
  - `rjsf-overrides.$defs`: the `type`/`properties` injection for all seven
    wrapped fields (this is the part that actually made wrapping work).
  - `rjsf-overrides.uiSchema.Home Range Method.home_range_args.args.anyOf`:
    label-hiding on each wrapper, matching the "inner field carries the
    visible title" rule from §3.
  - `param.yaml`, `test-cases.yaml`: field paths updated to the new nested
    shape (e.g. `speed_settings.max_speed_factor` instead of
    `max_speed_factor`).

## 7. Known remaining limitation

None for the field-blanking bug itself - the fix is complete and verified
(zero blank fields across a 7-step, double-loop, 3-variant cycle test
against the real compiled schema). Two unrelated, pre-existing cosmetic
items were surfaced during this work but are **not** fixable from schema
authoring alone (they'd require changes to ecoscope-desktop's own component
templates, out of scope here):

- Array/multi-select fields (e.g. "Percentile Levels") render with the same
  oversized heading style as objects - see §3's last note.
- The `ecoscope:advanced` flag does not appear to move a field into the
  "Advanced Configuration" accordion when the field lives inside an
  `anyOf`-selected branch's own properties (observed: `crs`/`percentiles`
  already carry this flag but still render as regular, always-visible
  fields at that nesting level).
