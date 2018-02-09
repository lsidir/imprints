library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

-- comparator component
--
--        in_config
--            |
--        /[ > ] \
--       /        +
--      / /[ > ] / \
--     / /          +  in_pos
--    / / /[ > ] \ / \   \
--   / / /        +   \   + out_pos
--  / / / /[ > ] /     \ /
-- in_data              +
--  \ \ \ \[ > ] \     /
--   \ \ \        +   /
--    \ \ \[ > ] / \ /
--     \ \          +
--      \ \[ > ] \ /
--       \        +
--        \[ > ] /
--

entity comparator is
generic (
	POS_WIDTH : integer := 6;     -- 6 bits for position 2^6 = 64 number of buckets
	DATA_WIDTH : integer := 64;   -- 64 bits data width
	LOG_NUM_COMP : integer := 3); -- log(number) of comparators
port (
	clk : in std_logic;           -- clokc
	resetn : in std_logic;        -- reset
	in_data : in std_logic_vector(DATA_WIDTH-1 downto 0);                      -- input data of DATA_WIDTH bits
	in_data_valid : in std_logic;                                              -- input data is valid
	in_pos : in std_logic_vector(POS_WIDTH-1 downto 0);                        -- input position to add
	in_config : in std_logic_vector(2**LOG_NUM_COMP*DATA_WIDTH-1 downto 0);    -- input configuration (bounds of bins in imprints)
	in_config_valid : in std_logic;                                            -- input configuration is valid

	out_data : out std_logic_vector(DATA_WIDTH-1 downto 0);                    -- output data is same as input
	out_data_valid : out std_logic;                                            -- output data valid same as input
	out_pos : out std_logic_vector(POS_WIDTH-1 downto 0));                      -- output new position (input position + new result)
end comparator;

architecture behavioral of comparator is

constant TREE_DEPTH : integer := LOG_NUM_COMP;                                 -- tree for adders
constant NUM_COMP : integer := 2**LOG_NUM_COMP;                                -- number of comparators

type bounds_type is array (NUM_COMP-1 downto 0) of std_logic_vector(DATA_WIDTH-1 downto 0); -- array for bounds of bins
signal bounds : bounds_type;

signal cmp_res : std_logic_vector(NUM_COMP-1 downto 0);
constant ZEROS : std_logic_vector(POS_WIDTH-2 downto 0) := (others => '0');

type intermediate_vector_type is array (NUM_COMP-1 downto 0) of unsigned(POS_WIDTH-1 downto 0); -- NUM_COMP x TREE_DEPTH+1 matrix of adder results
type intermediate_result_type is array (TREE_DEPTH downto 0) of intermediate_vector_type;
signal intermediate_result : intermediate_result_type;

signal in_data_valid_1d : std_logic;
signal in_data_valid_delay : std_logic_vector(TREE_DEPTH downto 0);

signal in_data_1d : std_logic_vector(DATA_WIDTH-1 downto 0);
type in_data_delay_type is array (TREE_DEPTH downto 0) of std_logic_vector(DATA_WIDTH-1 downto 0);
signal in_data_delay : in_data_delay_type;

type in_pos_delay_type is array (TREE_DEPTH-1 downto 0) of std_logic_vector(POS_WIDTH-1 downto 0);
signal in_pos_delay : in_pos_delay_type;


begin

GenX: for k in 0 to NUM_COMP-1 generate
	intermediate_result(0)(k) <= unsigned(ZEROS & cmp_res(k)); --concatenation at least significant
end generate GenX;

out_data_valid <= in_data_valid_delay(TREE_DEPTH);
out_data <= in_data_delay(TREE_DEPTH);

process(clk)
begin
if clk'event and clk = '1' then

	in_data_valid_1d <= in_data_valid; -- delay because of comparators
	in_data_valid_delay(0) <= in_data_valid_1d; -- delay tree depth
	for d in 1 to TREE_DEPTH loop
		in_data_valid_delay(d) <= in_data_valid_delay(d-1);
	end loop;

	in_data_1d <= in_data;
	in_data_delay(0) <= in_data_1d;
	for d in 1 to TREE_DEPTH loop
		in_data_delay(d) <= in_data_delay(d-1);
	end loop;

	in_pos_delay(0) <= in_pos;
	for d in 1 to TREE_DEPTH-1 loop
		in_pos_delay(d) <= in_pos_delay(d-1);
	end loop;

	if resetn = '0' then
		bounds <= (others => (others => '0'));
	else
		if in_config_valid = '1' then
			for i in 0 to NUM_COMP-1 loop
				bounds(i) <= in_config(DATA_WIDTH*(i+1)-1 downto DATA_WIDTH*i);
			end loop;
		end if;

		if in_data_valid = '1' then
			for i in 0 to  NUM_COMP-1 loop
				if in_data > bounds(i) then
					cmp_res(i) <=  '1';
				else
					cmp_res(i) <=  '0';
				end if;
			end loop;
		end if;

		for d in 0 to TREE_DEPTH-1 loop
			for k in 0 to 2**(TREE_DEPTH-d-1)-1 loop
				intermediate_result(d+1)(k) <= intermediate_result(d)(2*k) + intermediate_result(d)(2*k+1);
			end loop;
		end loop;

		out_pos <= std_logic_vector(unsigned(in_pos_delay(TREE_DEPTH-1)) + intermediate_result(TREE_DEPTH)(0));

	end if;
end if;
end process;

end behavioral;
