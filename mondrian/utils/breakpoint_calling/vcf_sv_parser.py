import pandas as pd
import vcf


class SvVcfData(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.reader = self._get_reader(filepath)
        self.caller = self._get_caller()

    def _get_caller(self):
        header_infos = self.reader.infos.values()

        if any(["lumpy" in str(v).lower() for v in header_infos]):
            return 'lumpy'
        elif any(["gridss" in str(v).lower() for v in header_infos]):
            return 'gridss'
        elif any(["svaba" in str(v).lower() for v in header_infos]):
            return 'svaba'
        else:
            raise Exception('unknown caller')

    @staticmethod
    def _get_reader(vcf_file):
        return vcf.Reader(filename=vcf_file)

    def _parse_vcf(self):

        for record in self.reader:

            data = {
                'chrom': record.CHROM,
                'pos': record.POS,
                'ref': record.REF,
                'alt': record.ALT,
                'qual': record.QUAL,
                'id': record.ID,
                'filter': record.FILTER,
            }

            # if data['filter'] == []:
            #     data['filter'] = None
            #
            # if data['filter'] is not None:
            #     continue

            info = record.INFO

            for k, v in info.items():
                if isinstance(v, list):
                    v = ';'.join(map(str, v))
                data[k] = v

            for sample in record.samples:
                sample_name = sample.sample
                sample_data = sample.data
                for k, v in sample_data._asdict().items():
                    if isinstance(v, list):
                        v = ';'.join([str(val) for val in v])
                    k = '{}_{}'.format(sample_name, k)
                    data[k] = v

            yield data

    def _group_bnds(self, calls):
        bnds = {}

        for record in calls:
            if self.caller == 'lumpy' and record['SVTYPE'] == 'INV':
                strands = record['STRANDS'].split(';')

                record['STRANDS'] = strands[0]
                yield (record,)
                record['STRANDS'] = strands[1]
                yield (record,)
            elif record['SVTYPE'] == 'BND':
                if 'MATEID' not in record:
                    continue

                if record['MATEID'] in bnds:
                    yield (record, bnds[record['MATEID']])
                    bnds.pop(record['MATEID'])
                else:
                    bnds[record['id']] = record
            else:
                yield record,

        assert len(bnds) == 0

    def _get_mates(self, records):
        ends_with_val = {
            'lumpy': '_1',
            'svaba': ':1',
            'gridss': 'h'
        }

        if records[0]['id'].endswith(ends_with_val[self.caller]):
            mate1, mate2 = records
        else:
            mate2, mate1 = records

        return mate1, mate2

    @staticmethod
    def _get_strand_from_alt(alt):
        """
        If the nucleotide comes first, then it is a "+" facing break at
        that site (see element "W" of the VCF4.2 specs, pg 13), otherwise
        it is "-" at that site. Then, if the bracket is ], then the
        partner breakpoint is "+", otherwise it is left facing
        :param alt:
        :type alt:
        :return:
        :rtype:
        """
        return alt[0].orientation

    def _get_strands(self, mate1, mate2):
        if self.caller == 'lumpy':
            strands_1 = mate1['STRANDS'].split(':')[0]
            strands_2 = mate2['STRANDS'].split(':')[0]
            assert strands_1 == strands_2[::-1]
            return strands_1[0], strands_2[0]
        else:
            strand_1 = self._get_strand_from_alt(mate1['alt'])
            strand_2 = self._get_strand_from_alt(mate2['alt'])

        return strand_1, strand_2

    @staticmethod
    def _process_lumpy_unmatched_record(record):
        record = record[0]

        strands = record['STRANDS'].split(':')[0]
        assert len(strands) == 2

        outdata = {
            'chromosome_1': record['chrom'],
            'position_1': record['pos'],
            'chromosome_2': record['chrom'],
            'position_2': record['END'],
            'strand_1': strands[0],
            'strand_2': strands[1],
            'type': record['SVTYPE']
        }

        return outdata

    def _process_bnd_call(self, record):

        assert len(record) == 2

        mate1, mate2 = self._get_mates(record)
        strand_1, strand_2 = self._get_strands(mate1, mate2)

        assert mate1['SVTYPE'] == mate2['SVTYPE']

        outdata = {
            'chromosome_1': mate1['chrom'],
            'position_1': mate1['pos'],
            'strand_1': strand_1,
            'chromosome_2': mate2['chrom'],
            'position_2': mate2['pos'],
            'strand_2': strand_2,
            'type': mate1['SVTYPE']
        }

        return outdata

    def _filter_low_qual_calls(self, calls):

        for call in calls:

            if len(call) == 1 and self.caller == 'lumpy':

                if call[0]['filter'] and 'LOW_QUAL' in call[0]['filter']:
                    continue
            else:
                assert len(call) == 2

                if call[0]['filter'] and 'LOW_QUAL' in call[0]['filter'] and 'LOW_QUAL' in call[1]['filter']:
                    continue

            yield call

    def fetch(self):
        records = self._parse_vcf()
        records = self._group_bnds(records)
        records = self._filter_low_qual_calls(records)

        for record in records:
            if len(record) == 1 and self.caller == 'lumpy':
                yield self._process_lumpy_unmatched_record(record)
            else:
                yield self._process_bnd_call(record)

    def as_data_frame(self):
        data = [record for record in self.fetch()]

        data = pd.DataFrame(data)
        data['caller'] = self.caller

        data['breakpoint_id'] = data.index

        data['breakpoint_id'] = data['breakpoint_id'].astype(str) + '_' + data['caller']

        return data
