#!/usr/bin/env bash
BIN=bin
PROG=../fmsi2
TESTS=testfiles

mkdir -p $BIN

# The commands of FMSI 1, on its index of several files.
LEGACY="$PROG legacy"

$LEGACY index $TESTS/integration_a.fa
$LEGACY index $TESTS/integration_b.fa

$LEGACY merge -p $TESTS/integration_a.fa -p $TESTS/integration_b.fa -r $BIN/merged.fa

$LEGACY query -k 3 -q $TESTS/queries.txt $TESTS/integration_a.fa > $BIN/a.txt 2> /dev/null
$LEGACY lookup -k 3 -q $TESTS/queries.txt $TESTS/integration_a.fa > $BIN/a_hash.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt -f xor $TESTS/integration_a.fa > $BIN/a_xor.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt $TESTS/integration_b.fa > $BIN/b.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt -f xor $TESTS/integration_b.fa > $BIN/b_xor.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt $BIN/merged.fa > $BIN/merged.txt 2> /dev/null
$LEGACY query -k 3 -q $TESTS/queries.txt -f xor $BIN/merged.fa > $BIN/merged_xor.txt 2> /dev/null

#$LEGACY normalize -k 3 -p $BIN/merged.fa -s -l > $BIN/merged_normalized.fa 2> /dev/null
$LEGACY normalize -k 3 -s $BIN/merged.fa > $BIN/merged_normalized2.fa 2> /dev/null



diff $TESTS/result_a_complements.txt $BIN/a.txt || exit 1
echo "a.txt OK"
diff $TESTS/result_a_complements_hash.txt $BIN/a_hash.txt || exit 1
echo "a_hash.txt OK"
diff $TESTS/result_b_complements.txt $BIN/b.txt || exit 1
echo "b.txt OK"
diff $TESTS/result_b_complements_xor.txt $BIN/b_xor.txt || exit 1
echo "b_xor.txt OK"
diff $TESTS/result_merged_complements.txt $BIN/merged.txt || exit 1
echo "merged.txt OK"

#diff $TESTS/result_normalized.txt $BIN/merged_normalized.fa
diff $TESTS/result_normalized2.txt $BIN/merged_normalized2.fa || exit 1
echo "merged_normalized2.fa OK"

# FMSI 2: one index file, with the same answers as FMSI 1 in every layout.
$PROG index -k 3 -o $BIN/a2 $TESTS/integration_a.fa 2> /dev/null || exit 1
[ -f $BIN/a2.fmsi ] || exit 1
$PROG info --check $BIN/a2.fmsi > $BIN/a2.info || exit 1
grep -q "^k	3$" $BIN/a2.info || exit 1
for layout in "" "-m plain" "-m packed" "-m rrr"; do
  for streaming in "" "-S"; do
    $PROG query $layout $streaming -q $TESTS/queries.txt $BIN/a2.fmsi > $BIN/a2.txt 2> /dev/null || exit 1
    diff $TESTS/result_a_complements.txt $BIN/a2.txt || exit 1
  done
done
# The index is also found by its name without the extension.
$PROG query -q $TESTS/queries.txt $BIN/a2 > $BIN/a2.txt 2> /dev/null || exit 1
diff $TESTS/result_a_complements.txt $BIN/a2.txt || exit 1
echo "FMSI 2 membership OK"

# Dictionary queries: the hashes are those of FMSI 1.
$PROG lookup -q $TESTS/queries.txt $BIN/a2.fmsi > $BIN/a2_hash.txt 2> /dev/null || exit 1
diff $TESTS/result_a_complements_hash.txt $BIN/a2_hash.txt || exit 1
# The type of a min-one mask is found out and recorded.
MINONE=data/GCF_test2_ones.fna.ms.k31
$PROG index -o $BIN/minone $MINONE 2> /dev/null || exit 1
$PROG info $BIN/minone.fmsi > $BIN/minone.info || exit 1
grep -q "^mask_type	min-one" $BIN/minone.info || exit 1
echo "FMSI 2 dictionary OK"

echo "All tests passed"

