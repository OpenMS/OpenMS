# OpenDIA issue map

[Master plan](README.md) · [Shared architecture](architecture.md)

Epic: [OpenMS/OpenMS#10661](https://github.com/OpenMS/OpenMS/issues/10661). Existing WP IDs are retained under five tracks. WP-09 was transferred from OpenMS-TimSim-Fixture#11 to OpenMS#10675, preserving its epic relationship.

| WP | Issue | Completion blockers |
|---|---|---|
| 00 Baseline | [#10662](https://github.com/OpenMS/OpenMS/issues/10662) | None |
| 01 Contracts | [#10663](https://github.com/OpenMS/OpenMS/issues/10663) | None |
| 02 Memory budget | [#10664](https://github.com/OpenMS/OpenMS/issues/10664) | 01 |
| 03 Raw tiles | [#10665](https://github.com/OpenMS/OpenMS/issues/10665) | 01 |
| 04 Library store | [#10666](https://github.com/OpenMS/OpenMS/issues/10666) | 01 |
| 05 Extraction/calibration adapter | [#10667](https://github.com/OpenMS/OpenMS/issues/10667) | 01 |
| 06 Retrieval | [#10668](https://github.com/OpenMS/OpenMS/issues/10668) | 03, 04, 05 |
| 07 Family scoring | [#10669](https://github.com/OpenMS/OpenMS/issues/10669) | 05, 06 for final scoring, none for early null/decoy experiment |
| 08 Output/FDR | [#10670](https://github.com/OpenMS/OpenMS/issues/10670) | 01 |
| 09 Integration | [#10675](https://github.com/OpenMS/OpenMS/issues/10675) | 02 through 08 for completion, starts on 01 |
| 10 Release gates | [#10671](https://github.com/OpenMS/OpenMS/issues/10671) | 00, 09 |

## Keeping docs and issues aligned

- Use these completion blockers for GitHub's native **blocked by** relationships. They do not prevent early experiments or adapter work. WP-01 closes on v0 contracts, not scientific contract freeze, so WP-07's review creates no dependency cycle.
- Keep the master short. Each WP has objective, scope, dependencies and 3–6 local merge checks. Shared policies live in architecture, whole-workflow gates in WP-09/10.
- Update existing issue bodies, not new duplicate issues. Preserve discussion, assignment and unchecked implementation status. Published issue links must be absolute issue URLs or pinned blob URLs, never pasted repository-relative links.
- Publish documentation to `singjc/OpenMS:rewrite/opendia` before linking its new revision. Use the resulting commit SHA for issue documentation links, so readers see exactly the reviewed version.
- Close packages with commands, pinned inputs/build and observed test results for their local checks. Keep end-to-end scientific/resource acceptance open in WP-09/10 until measured.
