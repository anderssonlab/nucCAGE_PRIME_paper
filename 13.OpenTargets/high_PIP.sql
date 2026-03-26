SELECT
  p.name AS prime_locus,
  p.facets AS prime_facets,
  p.score AS prime_score,
  cs.variantId,
  cs.pValueExponent,
  cs.beta,
  s.studyId,
  s.traitFromSource,
  t.approvedSymbol AS gene,
  locus_item.element.posteriorProbability AS pip_score,
  l2g.score AS l2g_score
FROM
  `andersson-lab.PRIME_FANTOM5.agnostic_distal_075` AS p
JOIN
  `open-targets-prod.platform.credible_set` AS cs
ON
  p.chromosome = cs.chromosome
  AND cs.position BETWEEN p.start AND p.end
  AND cs.studyType = "gwas"
JOIN
  `open-targets-prod.platform.l2g_prediction` AS l2g
ON
  cs.studyLocusId = l2g.studyLocusId
JOIN
  `open-targets-prod.platform.target` AS t ON l2g.geneId = t.id
JOIN
  `open-targets-prod.platform.study` AS s ON cs.studyId = s.studyId
CROSS JOIN UNNEST(cs.locus.list) AS locus_item
WHERE
  locus_item.element.variantId = cs.variantId
  AND locus_item.element.posteriorProbability > 0.3
ORDER BY 
  p.score * pip_score DESC