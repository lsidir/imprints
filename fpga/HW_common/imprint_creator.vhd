library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

entity imprint_creator is
generic(
	ADDR_WIDTH : integer := 32;
	DATA_WIDTH : integer := 64;
	COMP_WIDTH : integer := 8;
	POS_WIDTH : integer := 6;
	LOG_NUM_COMP : integer := 3);	-- 64 bits data width
port(
	clk : in std_logic;
	resetn : in std_logic;

	read_request : out std_logic;
	read_request_address : out std_logic_vector(ADDR_WIDTH-1 downto 0);
	read_request_almostfull : in std_logic;

	read_response : in std_logic;
	read_response_data : in std_logic_vector(511 downto 0);
	read_response_address : in std_logic_vector(ADDR_WIDTH-1 downto 0);

	write_request : out std_logic;
	write_request_address : out std_logic_vector(ADDR_WIDTH-1 downto 0);
	write_request_data : out std_logic_vector(511 downto 0);
	write_request_almostfull : in std_logic;

	write_response : in std_logic;

	start : in std_logic;
	done : out std_logic;

	-- Parameters
	number_of_CL_to_request : in std_logic_vector(31 downto 0));
end imprint_creator;

architecture structure of imprint_creator is

constant LOG2_REORDER_DEPTH : integer := 9;
constant MAX_ON_THE_FLY : integer := 2**(LOG2_REORDER_DEPTH-1);

signal NumberOfCacheLinesToRequest : unsigned(31 downto 0) := (others => '0');
signal NumberOfRequestedReads : unsigned(ADDR_WIDTH-1 downto 0) := (others => '0');
signal NumberOfCompletedReads : unsigned(31 downto 0) := (others => '0');
signal NumberOfRequestedWrites : unsigned(ADDR_WIDTH-1 downto 0) := (others => '0');
signal NumberOfCompletedWrites : unsigned(31 downto 0) := (others => '0');

signal reorder_start_address_adjust : std_logic;
signal reorder_start_address : std_logic_vector(ADDR_WIDTH-1 downto 0);
signal reordered_response_data : std_logic_vector(511 downto 0);
signal reordered_resonse : std_logic;

component reorder
generic(
	ADDRESS_WIDTH : integer := 32;
	LOG2_REORDER_DEPTH : integer := 8);
port (
	clk : in std_logic;
	resetn : in std_logic;

	start_address_adjust : std_logic;
	start_address : in std_logic_vector(ADDRESS_WIDTH-1 downto 0);
	in_trigger : in std_logic;
	in_address : in std_logic_vector(ADDRESS_WIDTH-1 downto 0);
	in_data : in std_logic_vector(511 downto 0);
	out_data : out std_logic_vector(511 downto 0);
	out_valid : out std_logic);
end component;

begin

reordering: reorder
generic map (
	ADDRESS_WIDTH => ADDR_WIDTH,
	LOG2_REORDER_DEPTH => LOG2_REORDER_DEPTH)
port map (
	clk => clk,
	resetn => resetn,

	start_address_adjust => reorder_start_address_adjust,
	start_address => reorder_start_address,
	in_trigger => read_response,
	in_address => read_response_address,
	in_data => read_response_data,
	out_data => reordered_response_data,
	out_valid => reordered_resonse);

process(clk)
begin

if clk'event and clk = '1' then
	NumberOfCacheLinesToRequest <= unsigned(number_of_CL_to_request);

	if resetn = '0' then
		NumberOfRequestedReads <= (others => '0');
		NumberOfCompletedReads <= (others => '0');
		NumberOfRequestedWrites <= (others => '0');
		NumberOfCompletedWrites <= (others => '0');

		reorder_start_address_adjust <= '0';
		reorder_start_address <= (others => '0');

		read_request <= '0';
		read_request_address <= (others => '0');

		write_request <= '0';
		write_request_address <= (others => '0');

		done <= '0';
	else -- resetn != '0'

		read_request <= '0';
		reorder_start_address_adjust <= '0';
		if start = '1' and NumberOfRequestedReads < NumberOfCacheLinesToRequest and read_request_almostfull = '0' and write_request_almostfull = '0' and NumberOfRequestedReads - NumberOfCompletedReads < MAX_ON_THE_FLY then
			read_request <= '1';
			read_request_address <= std_logic_vector(NumberOfRequestedReads);
			if NumberOfRequestedReads = 0 then
				reorder_start_address_adjust <= '1';
				reorder_start_address <= std_logic_vector(NumberOfRequestedReads);
			end if;
			NumberOfRequestedReads <= NumberOfRequestedReads + 1;
		end if;

		write_request <= '0';
		if reordered_resonse = '1' then
			NumberOfCompletedReads <= NumberOfCompletedReads + 1;

			write_request <= '1';
			write_request_address <= std_logic_vector(NumberOfRequestedWrites);
			write_request_data <= reordered_response_data;
			NumberOfRequestedWrites <= NumberOfRequestedWrites + 1;
		end if;

		if write_response = '1' then
			NumberOfCompletedWrites <= NumberOfCompletedWrites + 1;
		end if;

		if NumberOfCompletedWrites = NumberOfCacheLinesToRequest and NumberOfCacheLinesToRequest > 0 then
			done <= '1';
		end if;

	end if; -- resetn = '0'
end if; -- clk'event and clk = '1'
end process;
end structure;
